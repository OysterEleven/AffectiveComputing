function GNE = gne(x,Fs,F0,Tind)
%
% Finds the glottal-to-noise excitation ratio (GNE) of a speech signal
%
% Usage: ratio = gne(x,Fs)
%        ratio = gne(x,Fs,F0,Tind)
%
%      x: input signal
%     Fs: sampling frequency
%     F0: 1xN vector contaning F0 contour (optional)
%   Tind: Nx2 vector containing indices of analysis frames (optional)
%         If F0 or Tind are missing, gne calls getF0 to compute them.
%         If they're already available from a previous computation, pass
%         to speed up computation.
%
%  ratio: a 1xN vector containing the glottal-to-noise excitation ratio
%         for voiced frames. Voiceless frames are encoded as zeros.
%         This is based on the algorithm described in: 
%               Michaelis, D., Gramss, T., and Strube, H.W., 
%               "Glottal-to-Noise Excitation Ratio -- A New Measure for 
%               Describing Pathological Voices." Acustica. Vol. 83
%               pp. 700-706, 1997.
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



if nargin < 3
  vflag = 0;
else
  vflag = 1;
end  
  

% Downsample to 8kHz
if Fs > 8000
  Fsd = 8000;
  ratio = Fsd/Fs;
  s = resample(x,round(ratio*1000),1000); 
  if vflag; Tind = ceil(Tind*Fsd/Fs); end
else
  s = x;
  Fsd = 8000;
end

% Inverse filter to obtain residual and work with residual instead
win_length = round(30e-3*Fsd);		%Analysis window of 30 msecs.
win_step   = round(10e-3*Fsd);		%Shift of 10 msecs. between frames
P = 10;					%LPC order
warning off;
[ar,e,k] = lpcauto(s,P,[win_step win_length]);
warning on;
are = ar.*(sqrt(e)*ones(1,P+1));
u = lpcifilt(s,ar,k(:,1));

%%
%% If no voicing info passed, find voiced segments
%% and restrict analysis to the voiced frames from the
%% output of the pitch tracker. Otherwise, limit the
%% analysis to the samples corresponding to the analysis
%% indices passed
%%

if (~vflag)
  % Get pitch contour, strenghts, etc.
  [F0,strength,Tind] = getF0(s,Fsd); 
end
vind = find(F0 ~= 0);
numframes = length(F0);
GNE = zeros(1,numframes);



for n=1:length(vind)
  framenum = vind(n);
  to = Tind(framenum,1);
  tn = Tind(framenum,2);
  analysis_region = u(to:tn); 

  %% Find Hilbert envelopes from the frequency domain
  % Find spectrum of voiced excitations
  fftsize = 2^nextpow2(length(analysis_region));
  analysis_region = analysis_region(:)';
  U = fft(analysis_region,fftsize);

  % Band-filter using Hanning windows for bands 0-2kHz; 1-3kHz; and 2-4kHz
  % Only positive frequencies considered. Find Hilbert envelope by filtering in
  % frequency --> taking inverse FFT --> finding
  % absolute value of complex-valued signal
  Fnyq = Fsd/2;
  Bnyq = fftsize/2;
  mult = Bnyq/Fnyq;
  w1 = [1:2000*mult];
  w2 = [1000*mult:3000*mult];
  w3 = [2000*mult:4000*mult];
  U1 = zeros(1,fftsize); U2=U1; U3=U1;
  U1(w1) = U(w1).*hanning(length(w1))';
  U2(w2) = U(w2).*hanning(length(w2))';
  U3(w3) = U(w3).*hanning(length(w3))';
  u1 = abs(ifft(U1)); u2=abs(ifft(U2)); u3=abs(ifft(U3));

  % Determine pairwise cross-correlations between each of the
  % Hilbert envelopes for lags between -0.3ms<lag<0.3ms.
  % Normalize the xcorrs by the sqrt of the product of
  % the energy of each signal to obtain a number 0<n<1.
  % The glottal to noise excitation ratio is the max of
  % the max of each cross-correlation
  maxlag = ceil(0.3e-3*Fsd);
  x12 = xcorr(u1,u2,maxlag)/(sqrt(sum(u1.^2)*sum(u2.^2)));
  x13 = xcorr(u1,u3,maxlag)/(sqrt(sum(u1.^2)*sum(u3.^2)));
  x23 = xcorr(u2,u3,maxlag)/(sqrt(sum(u2.^2)*sum(u3.^2)));
  %figure; 
  %subplot(311); plot(x12)
  %subplot(312); plot(x13)
  %subplot(313); plot(x23)
  GNE(framenum) = max(max([x12' x13' x23']));
end


