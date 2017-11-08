function [LFtheta,LFsignal,LFerror,ghat,gtap] = glottal_invfilt(x,Fs,landpts)
%
% [LFtheta,LFsignal,LFerror,ghat,gwin] 
%                              = glottal_invfilt(x,Fs, glottal_landmarks)
%
% High-level function to identify phonatory cycles in a speech signal, 
% inverse-filter and fit the parametric LF model to each cycle.
%
% Inputs:
%	           x:	speech signal
%                 Fs:	sampling frequency
%  glottal_landmarks:   Nx4 array containing glottal phase details. The output
%                       of glottmarks.m (see GLOTTMARKS). Pass if already
%                       available to speed up computation, or default to empty.
%
% Outputs:
%
%            LFtheta:   a 4xN matrix, where each column contains the LF 
%                       parameters of a phonatory cycle
%           LFsignal:   a synthesized LF signal from LFtheta (same length as x)
%            LFerror:   a 2xN matrix with each column containing the MS-error 
%                       of the LF fit throughout the open and closed phases 
%                       respectively
%   		ghat:   the inverse filtered volume velocity waveform;
%               gwin:   a post-windowed version of ghat, with the closed 
%                       phase further de-emphasized.
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
 	


P = 18;	% (max.) order of LPC model to fit for inverse filtering


if nargin <3
  % Parse the speech waveform into "landmark" points
  landpts = glottmarks(x,Fs);
end
[N,dim2] = size(landpts);
CP_delta = landpts(:,4)-landpts(:,3);



%% For each glottal cycle
% - fit an AR model (using DAP) over the closed phase
% - fix the polynomial to get rid of "extraneous" poles
% - inverse filter the entire cycle
% - fit an LF model to the glottal estimate
Fr = [1/4 1/2 3/4 1];
Wts = [5 1  1 1];
yx = 0.1;
y0 = 1;

ghat = zeros(1,length(x));
gtap = zeros(1,length(x));
gLF  = zeros(1,length(x));
gLF2 = zeros(1,length(x));
cum_win = zeros(1,length(x));

% Process blocks of three cycles: the current one, the previous
% and the next. Retain only the results for the current cycle, and
% advance the counter by one.
pcd = 0.9;
disp('Inverse filtering and fitting an LF model ....')
disp('to each phonatory cycle ......................')
disp('(This takes time. Go get a cup of coffee).....')
for n=1:N
   if n==round(0.25*N);
     disp('25% completed ................................');
   elseif n==round(0.5*N)
     disp('50% completed ................................');
   elseif n==round(0.75*N)
     disp('75% completed ................................');
   end
   % Segment three adjacent cycles, and build vector with closed-phase
   % indices. Treat the boundary cases (n=1,N) with only two cycles.
   if (n==1)
     segn = x(landpts(n,1):landpts(n+1,4));
     offset = landpts(n,1)-1;
     cp1=[];
     cp2=[landpts(n,3):landpts(n,3)+round(pcd*CP_delta(n))]-offset;
     cp3=[landpts(n+1,3):landpts(n+1,3)+round(pcd*CP_delta(n+1))]-offset;
   elseif (n==N)
     segn = x(landpts(n-1,1):landpts(n,4));
     offset = landpts(n-1,1)-1;
     cp1=[landpts(n-1,3):landpts(n-1,3)+round(pcd*CP_delta(n-1))]-offset;
     cp2=[landpts(n,3):landpts(n,3)+round(pcd*CP_delta(n))]-offset;
     cp3=[];
   else
     segn = x(landpts(n-1,1):landpts(n+1,4));
     offset = landpts(n-1,1)-1;
     cp1=[landpts(n-1,3):landpts(n-1,3)+round(pcd*CP_delta(n-1))]-offset;
     cp2=[landpts(n,3):landpts(n,3)+round(pcd*CP_delta(n))]-offset;
     cp3=[landpts(n+1,3):landpts(n+1,3)+round(pcd*CP_delta(n+1))]-offset; 
   end
   cp = [cp1 cp2 cp3];	%Closed-phase indices

   % Adaptive pre-emphasis
   alpha = real(lpc(segn,1));
   seg_emph = filter(alpha,1,segn);

   % Extract the closed phase portion from the pre-emphasized signal
   closed_phase_seg = [seg_emph(cp1)' seg_emph(cp2)' seg_emph(cp3)'];

   % Apply DAP modeling to the data from the 3 adjacent close phases
   ARcf{1,n} = dapw(closed_phase_seg,P,100,200,Fr,Wts);
   % If numerical problems arise, decrease P until it's stable
   unstable = sum(isnan(ARcf{1,n})>0) | sum(isinf(ARcf{1,n})>0);
   while (unstable & (P>= 2))
      str=['Unstable DAP fit. Decreasing order from ' num2str(P) ' to ' num2str(P-1)];
      disp(str);
      P = P-1;
      ARcf{1,n} = dapw(closed_phase_seg,P,100,200,Fr,Wts);
      unstable = sum(isnan(ARcf{1,n})>0) | sum(isinf(ARcf{1,n})>0);
      if ~unstable disp('Stable fit found.'); end
   end

   % Construct inverse filter polynomial by "fixing" the roots
   if ~unstable
     ARif{1,n} = formant_reshape(ARcf{1,n},Fs);
   else
     disp('No stable fit found. Proceeding')
     ARif{1,n} = ones(1,P);
   end

   % Inverse-filter to filter out effect of closed phase
   g = filter(ARif{1,n},1,segn);

   %%% Post-multiply after filtering
   ph1 = landpts(n,2)-landpts(n,1)+1;
   ph2 = landpts(n,4)-landpts(n,2);

   % Window segment for the closed phase (decaying exponential)
   a = -log(yx/.5)/ph2; %Calculate rate of decay of exponential
   wcp = .5*exp(-a*[0:1:ph2-1]);
   
   % Window segment for the open phase (logistic fxn. with sharp rise)
   xm = 0.1*ph1; b = log((1-yx)/yx)/xm;
   wop = 1./(1+exp(-b*([0:1:ph1-1]-xm)));

   % String the two portions together
   postwin = [wop wcp];	%Post-processing ad-hoc window
   inds = [landpts(n,1):landpts(n,4)]-offset;
   g_center = g(inds);
   g_c_win = g_center.*postwin';

   ghat(landpts(n,1):landpts(n,4))=g_center';
   gtap(landpts(n,1):landpts(n,4))=g_c_win';
   cum_win(landpts(n,1):landpts(n,4))=postwin';

   %[th,gth] = fitLF(g_c_win);
   %gth_cum{1,n}=gth; th_cum(:,n)=th';
   %gLF(landpts(n,1):landpts(n,4))=gth';

   Neval = landpts(n,2)-landpts(n,1)+1;
   %[th2,gth2] = fitLF(g_center,Neval,'off');
   [th2,gth2,eop,ecp] = fitLF(g_center,[],'off');
   gth2_cum{1,n}=gth2; th2_cum(:,n)=th2'; error_cum(:,n) = [eop ecp]';
   gLF2(landpts(n,1):landpts(n,4))=gth2';
end

LFtheta = th2_cum;
LFerror = error_cum;
LFsignal = gLF2;












%t = [1:length(x)]/Fs;
%figure
%subplot(411)
%plot(t,x)
%%axis([.77 .81 -1 1])
%axis([1.305 1.335 -.5 .5])
%title('(a)')
%subplot(412)
%plot(t,ghat,'r')
%%axis([.77 .81 -1.5 1])
%axis([1.305 1.335 -1 1])
%title('(b)')
%subplot(413)
%plot(t,cum_win,'r--')
%hold on
%plot(t,gtap,'k')
%%plot(t,gLF)
%%%axis([.77 .81 -1.5 1])
%axis([1.305 1.335 -1 1])
%%title('(c)')
%subplot(414)
%plot(t,gLF)
%hold on
%plot(t,gLF2,'r--')
%%axis([.77 .81 -1.5 0.5])
%axis([1.305 1.335 -1 1])
%xlabel('Time (secs.)')
%title('(d)')
%
%
