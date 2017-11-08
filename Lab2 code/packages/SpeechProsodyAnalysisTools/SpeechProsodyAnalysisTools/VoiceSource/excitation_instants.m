function [instants,phaseslope,res,signal]=excitation_instants(signal,original_freq,avg_pitch_freq,downsampling_freq, multiplier,smoothwidth,LPCorder,LPCframedetails,verbose)
%%%%%%%%%%%%%%%%%% EXCITATION_INSTANTS.M %%%%%%%%%%%%%%%
%
% Determine the instants of significant excitation (e.g.
% glottal closures) in a speech signal, using the positive
% zero-crossings of the phase slope function of the signal's
% linear-prediction residual.
%
%                  ------ usage ------
%
%[instants,phaseslope,res,signal]=excitation_instants(signal,original_freq,
% +avg_pitch_freq+,+downsampling_freq+,+multiplier+,+smoothwidth+,+LPCorder+,
% +LPCframedetails+);
% 
% where +...+ indicates an optional (defaulted) argument.
%
%
%                 ------ inputs ------
% signal: Input signal for which excitation instants should be found. 
% original_freq: Sampling frequency of input signal. 
% avg_pitch_freq: Assumed approximate fundamental pitch frequency,
%          in Hz. Obviously this should be higher for female speakers
%          than for male speakers.  (default: 150)
% downsampling_freq: downsample audio signal to this frequency 
%          before proceeding (default: 8000 Hz)
% multiplier: Specifies the length of the short-time Fourier transform
%          analysis window in multiples of the fundamental pitch period.
%          This should be between 1 and 2. (default: 1)
% smoothwidth: Width of the Hanning window used to smooth the pre-
%          liminary phaseslope plot. (default: 30)
% LPCorder:Order of LPC analysis to be calculated to determine the
%          residual. (default: 10)
% LPCframedetails: [framespacing framewidth skip]
%          Perform aforementioned LPC analysis using frames of  
%          width framewidth, centered framespacing samples apart,
%          skipping skip samples at the beginning.   
%          (default: [80 200 0])
% 
% 
%                 ------ outputs ------
%
% instants: Vector, variable length. Indices of all instants of
%          significant excitation as determined by the algorithm.
% phaseslope: Vector, same length as signal.  The nth entry is 
%          the slope of the best linear fit to the unwrapped 
%          phase of a short-time Fourier transform centered at
%          point n in the original signal. 
%     res: Residual of LPC analysis. 
%  signal: Original signal for which excitation instants were found.
%
% 
% Important notes: This updated implementation performs the excitation
%          instants calculation on a downsampled version of the signal
%          that has Fs=downsampling_freq, but then expands the resulting
%          vector of excitation instants to correspond to the original 
%          signal. If you manually expand the output it will have been
%          done twice.
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



%%%%%% input handling %%%%%%%%%%%%%%%%%%%%%%%%%%'
% 
%
if(nargin<9)
  verbose = 0;
end
if(nargin<8)
     LPCframedetails=[80 200 0];
end
if(nargin<7)
     LPCorder=10;
end
if(nargin<6)
     smoothwidth=30;
end
if(nargin<5)
     multiplier=1;
end
if(nargin<4)
     downsampling_freq=8000;
end
if(nargin<3)
     avg_pitch_freq=150;
end
if(nargin<2)
     error('The input signal AND its original sampling frequency are required.') 
end
%
%
if(size(original_freq)~=[1 1]) error('original_freq must be a scalar.'); 
end
if(size(avg_pitch_freq)~=[1 1]) error('avg_pitch_freq must be a scalar.'); 
end
if(size(downsampling_freq)~=[1 1]) error('downsampling_freq must be a scalar.');
end
if(size(multiplier)~=[1 1]) error('multiplier must be a scalar.'); 
end
if(size(smoothwidth)~=[1 1]) error('smoothwidth must be a scalar.'); 
end
if(size(LPCorder)~=[1 1]) error('LPCorder must be a scalar.'); 
end
if(size(LPCframedetails)~=[1 3]) error('LPCframedetails must be 1x3.'); 
end
%
%
%%%%% lpc analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
downsampled=downsample_to_n_hz(signal,original_freq,downsampling_freq);
if verbose
   fprintf(['\nEXCITATION_INSTANTS> Original signal'  ...
       ' downsampled from ' num2str(original_freq) ' to ' ...
        num2str(downsampling_freq) ' Hz.\n']);
end
tic;
[ar,e,k,res]=lpc_residual(downsampled,LPCorder,LPCframedetails,'hanning',0,0);
residual_time=toc;
if verbose
  fprintf(['EXCITATION_INSTANTS> Order ' int2str(LPCorder) ' LPC residual calculated in ' num2str(residual_time) ' sec.\n']);
end
%
%%%%% short-time fourier transform %%%%%%%%%%%
%
%  preliminary calculations
%
avg_pitch_per=1/avg_pitch_freq; % sec
if verbose
  fprintf(['EXCITATION_INSTANTS> Assumed pitch frequency is ' num2str(avg_pitch_freq) ' Hz; average period ' num2str(1000*avg_pitch_per) ' msec.\n']);
end
frametime=multiplier*avg_pitch_per; % sec
framesamples=floor(frametime/(1/downsampling_freq));
if (mod(framesamples,2)==0)
     framesamples=framesamples-1;
end%if
windowlength=framesamples;
if verbose
  fprintf(['EXCITATION_INSTANTS> Multiplier is ' num2str(multiplier) ', so STFT window is ' num2str(frametime*1000) ' msec, or ' num2str(windowlength) ' samples.\n']);
end
hannwin=window('hanning',windowlength);
if verbose
  fprintf('EXCITATION_INSTANTS> Ready for STFT calculation. Rock it down.\n\n');
end
%
%  real deal
%
tic;xform=stft(res,hannwin);stft_time=toc;
if verbose
  fprintf(['\nEXCITATION_INSTANTS> STFT completed in ', num2str(stft_time), ' sec. Calculating phase angles...\n']);
end
%
%%%%% phase angles %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
tic;phi=angle(xform);phase_time=toc;
if verbose
  fprintf(['EXCITATION_INSTANTS> Phase angles found in ' num2str(phase_time) ' sec. Unwrapping phase...\n']);
end
%
%%%%% unwrapped phase %%%%%%%%%%%%%%%%%%%%%%%%
%
tic;unwrapped=unwrap(phi,[],1);unwrap_time=toc;
%
if verbose
  fprintf(['EXCITATION_INSTANTS> Unwrapping done in ' num2str(unwrap_time) ' sec. Ready for regression...\n']);
end
%
%%%%% linear regression %%%%%%%%%%%%%%%%%%%%%%
%
% ugly, but 12x faster than regress 
% (plus THIS actually works)
%
[Q,R]=qr([1:size(unwrapped,1)]',0);
Qtrans=Q';
phaseslope=ones(1,size(unwrapped,2));
%
% perform a linear fit for each column
%
tic
for i=1:size(unwrapped,2)
phaseslope(i)=R\(Qtrans*unwrapped(:,i));
end%for
regress_time=toc;
%
% zero-mean it
%
phaseslope=phaseslope-mean(phaseslope);
if verbose
  fprintf(['EXCITATION_INSTANTS> Regression completed in ' num2str(regress_time) ' sec.\n']);
  fprintf(['EXCITATION_INSTANTS> Smoothing phase slope with size ' num2str(smoothwidth) ' Hanning window.\n']);
  fprintf('EXCITATION_INSTANTS> Zero-centering phase slope and generating plot...\n');
end
%
%%%%% phaseslope smoothing %%%%%%%%%%%%%%%%%%%%%
%
% use Hanning window smoothing
%
phaseslope=filter(hanning(smoothwidth),1,phaseslope);
%
% correct for the observed delay
% this gets within one sample for practical smoothing widths
%
approxdelay=floor((smoothwidth-2)/2); 
phaseslope=phaseslope(approxdelay+1:end);
%
%%%%% zero crossings %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% find instants of positive zero-crossings
%
unexpanded_instants=find(diff(sign(phaseslope))==2);
instants=unexpanded_instants*round(original_freq/downsampling_freq);
excitationvector=zeros(1,size(unwrapped,2));
excitationvector(instants)=0.15;
%
%
%%%%% outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% plot 'em %'
%
%figure 
%hold on
%plot(signal)
%stem(excitationvector,'r')
%title(['Instants of significant excitation calculated at Fs=' num2str(downsampling_freq) ', re-expanded to ' num2str(original_freq)]);%
%xlabel('Samples');
%ylabel('Signal: blue, excitation instants: red');
%zoom on
%
fprintf(['EXCITATION_INSTANTS> ' int2str(length(instants)) ' instants of significant excitation.\n']);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









