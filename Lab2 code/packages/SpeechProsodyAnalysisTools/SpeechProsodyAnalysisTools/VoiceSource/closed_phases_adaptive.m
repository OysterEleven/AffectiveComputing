function [closed,formantvec]=closed_phases_adaptive(signal,ex_instants,order,Fs,k,lmult,rmult,maxfreq,Nskip,verbose,ploton);
% CLOSED_PHASES_ADAPTIVE: Find the closed glottal phases of a speech signal 
% using Quatieri/Plumpe tracking of first-formant modulation. 
%
% Usage:
%
% closed=closed_phases_adaptive(signal,ex_instants,+order+,+Fs+,+k+,+lmult+,
%                               +rmult+,+maxfreq+,+Nskip+,+verbose+,+ploton+); 
% 
% where +argument+ designates optional (defaulted) arguments. Optional arguments
% may be skipped and defaulted by assigning them to the empty matrix []
%     Ex: cp = closed_phases_adaptive(signal,excits,[],8000,[],[],[],[],[],1) 
%  
% signal:  	input speech signal
% ex_instants:	a vector of signal's instants of significant
%        	excitation, (HELP ESCITATION_INSTANTS)
% order:	order for the sliding LPC covariance used to find input's
% 		first formant frequency (unless this order would cause 
%       	matrix singularity in a particular analysis window, in which 
%       	case the smallest numerically sufficient order is
%       	substituted for that window only. default:16);
% Fs:		input's sampling rate in Hz (default:16000);
% k:		controls the seeding of closed phase regions. The flattest 
%       	length-k segment of the search region of a formant track
% 		is selected as the closed phase region to be grown (def:5);
% lmult 
%  and 
% rmult: 	define the number of standard deviations from the mean a
%		formant point must fall within to be included in the closed 
%		phase region during leftward and  rightward growth, 
%		respectively. (Default: lmult=2 and rmult=2)
% maxfreq:	maximum pitch frequency in Hz.  Any excitation spikes
%       	in ex_instants corresponding to a higher frequency will 
%		be  removed (Default:550);
% Nskip:	number of samples by which to advance the sliding LPC 
%		covariance analysis window for each frame (default:1);
% verbose:	a boolean argument (0=off; 1=true) to enable stats dump on 
%		the Matlab environment at run time (default:0);
% ploton: 	a boolean argument (0=off; 1=true) to enable plotting of 
%		some intermediate results (useful when debugging) (default:0)
%
% ...and
%
% closed is a matrix of indices describing the location of the excitation
% (some spurious ones may be filtered out while running this code) as well 
% as the  beginning and end of the closed phases corresponding to each given 
% excitation instant:
%
%     [excit1 beginningofclosed1 endofclosed1;
%      excit2 beginningofclosed2 endofclosed2;
%         ...         ...           ...
%      excitn beginningofclosedn endofclosedn]
%
% Two final notes: (1) The search area, which is currently set to 
% 1/4 the length of any formant track, may easily be set on line 
% 127, search_area.. (2) Whether or not to include the P samples at 
% the beginning of a window as part of the speech used may easily be 
% set on line 146 with the subtraction of orders_actually_used(i).  
%
%

%
% Copyright (C) 2006  Raul Fernandez | galt@media.mit.edu
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You can obtain a copy of the GNU General Public License from
% http://www.gnu.org/licenses/gpl.txt or by writing to
%   The Free Software Foundation, Inc.,
%   51 Franklin Street, Fifth Floor,
%   Boston, MA  02110-1301, USA.
%





%%%%% INPUT HANDLING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
if (nargin<2); error('Not enough arguments!'); end
if (nargin<3 | isempty(order)); order=16; end;
if (nargin<4 | isempty(Fs)); Fs=16000; end ;  
if (nargin<5 | isempty(k)); k=5; end;
if (nargin<6 | isempty(lmult)); lmult=2; end
if (nargin<7 | isempty(rmult)); rmult=2; end;
if (nargin<8 | isempty(maxfreq)); maxfreq=550; end
if (nargin<9 | isempty(Nskip)); Nskip=1; end
if (nargin<10 | isempty(verbose)); verbose=0; end
if (nargin<11 | isempty(ploton)); ploton=0; end

%
%
%%%%% DO IT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Preprocess ex_instants to filter out ridiculously small intervals
%
smallest_period=1/maxfreq; % seconds
smallest_window=floor(Fs*smallest_period); % samples
ex_instants_filt=[];
for i=1:length(ex_instants)-1
if ( (ex_instants(i+1)-ex_instants(i)) > smallest_window)
     ex_instants_filt=[ex_instants_filt ex_instants(i)];
end%if 
end%for
if verbose
  fprintf(['\nCLOSED_PHASES_ADAPTIVE> ' ... 
	num2str(length(ex_instants)-length(ex_instants_filt)) ...
        ' spurious excitation instants removed before processing.\n']);
end

%
% For every consecutive pair of excitation spikes:
%
if (length(ex_instants_filt) <2 )	%Degenerate case: there's only one
  closed = [];				%or no spikes
  formantvec = [];
  return;
end

for i=1:(length(ex_instants_filt)-1) 
this_spike=ex_instants_filt(i);  
next_spike=ex_instants_filt(i+1); 
distance=next_spike-this_spike; 
Nw=floor(distance/4);
%
% Reduce order if necessary to prevent LPC failure
%
orders_actually_used(i)=order; % try to use default order
if (orders_actually_used(i)>Nw-3) % possible Cholesky matrix singularity!
  orders_actually_used(i)=Nw-3; % so reduce to maximum allowable order
  if verbose
    fprintf(['CLOSED_PHASES_ADAPTIVE> Order ' num2str(order) ' was reduced to ' ...
					   num2str(orders_actually_used(i)) ...
					  ' between instants ' num2str(i) ...
					  ' and ' num2str(i+1) ' for Nw=' ...
					  num2str(Nw) '.\n']);
  end
end
%
% Get N-Nw+2 formant values
%

if orders_actually_used(i) == 1;   % An order of 1 does not yield any formants. Run with
  orders_actually_used(i)=2;       % order of at least 2.
end

[formants{i},Tcell{i}]=formant1a(signal(this_spike:next_spike-1), ... 
				Nw,Nskip,orders_actually_used(i),Fs); % gets N-Nw+2 formants
%
% Call upon find_flattest_segment to seed and grow closed phase regions
%
%search_area=distance-Nw-5;
search_area=floor(.25*length(formants{i}));
[first,last,seed]=find_flattest_segment((formants{i}(1:search_area)),k,lmult,rmult);
closed_frames(i,1)=first; % save segment beginnings
closed_frames(i,2)=last; % save segment endings 
seeds(i)=seed; % save seeds
end%for
%
%
%%%%% STABLE FRAMES -> STABLE SAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Above we calculated the LPC covariance FRAMES over which the 
% first formant is stable. However, to properly perform the inverse 
% filtering, we need to know the samples of the original speech 
% signal corresponding to these frames. This correspondence is kept
% by the T matrix output by formant1.m. So...
%
NumClosedFrames = size(closed_frames,1);	%Bug fixed on 6.23.03 -- raul
                                                %Code was using length, instead of size(.,1)
                                                %and erring out when closed frames was only one row 
closed=zeros(NumClosedFrames,2);
for i=1:NumClosedFrames
     % MATCHING WITH SIGNAL
     closed(i,1)=Tcell{i}(closed_frames(i,1),1)+ex_instants_filt(i)-orders_actually_used(i); 
                                                               % use those P samples for now 
     closed(i,2)=Tcell{i}(closed_frames(i,2),2)+ex_instants_filt(i); 
end

% Augment the closed phase matrix with the excitation instants that were
% not removed
closed = [ex_instants_filt(1:end-1)' closed];

%
%
%%%%% PLOT OUTPUT %%%%%%
%
%
if ploton
  outputvec=zeros(1,length(signal));
  instantsvec=zeros(1,length(signal));
  seedsvec=zeros(1,length(signal));
  formantvec=zeros(1,length(signal));
  for i=1:length(closed)
       formantvec(ex_instants_filt(i):(ex_instants_filt(i)+length(formants{i})-1))=...
                    (formants{i}/3000);
       seedsvec(seeds(i)+ex_instants_filt(i))=.20;
  end;
  instantsvec(ex_instants_filt)=.25;
  for i=1:length(closed)
    outputvec((closed_frames(i,1)+ex_instants_filt(i)):...
                        (closed_frames(i,2)+ex_instants_filt(i)))=.15;
  end;
  figure;
  plot(signal,'b');
  plot(formantvec,'m');
  hold on;
  stem(outputvec,'r');
  stem(seedsvec,'g');
  stem(instantsvec,'k');
  zoom on;     
  title('Closed Phase Frames (not samples)')
  xlabel('r:closed phase frames, g:closed phase seeds, bk:excitations, m:formant tracks');
end













