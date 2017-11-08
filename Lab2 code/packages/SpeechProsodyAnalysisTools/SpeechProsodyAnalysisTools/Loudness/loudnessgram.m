function [N,rmsvals,T_ind,Barkscale,N_spec] = loudnessgram(sp,Fs,A0,win_dur,time_step)
%
% [N,rmsvals,Tind,Barks,N_spec] = loudnessgram(sp,Fs,A0,win_dur,win_step)
%
% 	 sp:  speech waveform
%	 Fs:  sampling rate
%	 A0:  RMS level to normalize the speech waveform before calling
%	      loudness algorithm ~ (0.25-0.5). Default is 0.25.
%  wind_dur:  duration of analysis window (in secs.) Default is 40e-3
% time_step:  length of step (in secs.) between adjacent windows. 
%             (Default is 10e-3);
% 
% Outputs:
%  	  N:  loudness curve; each sample is the absolute loudness of a
%	      short-time segment of the normalized waveform, calculated
%	      using the DIN 45631 loudness standard
%      Tind:  time indices of each segment used in calculating N
%     Barks:  Bark scale in 0.5 bark increment
%    N_spec:  a length(Barks)xlength(Tind) matrix with each column 
%             representing the specific loudness (the distrib. of loudness 
%             over bark freq.) for each short-time segment. N(k) is the 
%             integral over barks of N_spec(:,k)
%
% Refs: -"Psychoacoustics," Zwicker, E. and Fastl, H. Springer Verlag
%       -"Program for Calculating Loudness According to DIN 45631 (ISO 532B),"
%        Zwicker, E., Fastl, H. et al. J. Acout. Soc. Jpn (E). 
%        Vol. 12 No. 1, 1991
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



if (nargin < 3 | isempty(A0)) A0 = 0.25; end 
if (nargin < 4 | isempty(win_dur)) win_dur = 40e-3; end
if (nargin < 5 | isempty(time_step)) time_step = 10e-3; end


% Check that Fs is a 'valid' sampling rate 
% (i.e., handled by the routine that generates filters)
Fsallowed = [4000 8000 16000 22050 44100];
if isempty(intersect(Fs,Fsallowed)); 
  error('Supported sampling rates are 4, 8, 16, 22.05 and 44.1kHz. Resample before calling this routine.')
end


% Fix the RMS value of entire signal to A0
sp = A0*sp/rms(sp);

% Parameters
win_length = round(win_dur*Fs);
win_step = round(time_step*Fs);
L = length(sp);
kmax = floor(1+(L-win_length)/win_step);
psdNoverlap = 0;


% Generate filters
Nfft = 512;
df = (Fs-0)/(Nfft);
freqs = [0:df:Fs/2];
[H, err]=myGenerateFilters(freqs,Fs);
H = 0.9639*H;

N = []; rmsvals=[]; N_spec=[];
% Bark ranges (to transform from 1/10 bark increments to half-bark increments)
Barks = [1:240]/10;
barkrange = [0:0.5:24];
Barkscale = [.25:.5:24];
Nbrks = length(barkrange)-1;
for k=1:Nbrks
  B1 = find(Barks >= barkrange(k)); b1(k) = B1(1);
  B2 = find(Barks < barkrange(k+1)); b2(k) = B2(end);
end

% Short-term analysis
for k=1:kmax
  % Choose a segment and store the index locations (onset and offset)
  seg = sp((k-1)*win_step+1:(k-1)*win_step+win_length);
  T_ind(k,:) = [(k-1)*win_step+1 (k-1)*win_step+win_length];

  % Find the segment's power spectral density
  Xseg = psd(seg,Nfft,Fs,Nfft,psdNoverlap);

  % Convert to dBs (taking care of numerical zeros)
  zeroind = find(Xseg == 0);
  Xseg(zeroind) = eps;
  alpha = 1;
  XdB = 10*log10(Xseg/(alpha^2));
  XdBn = XdB';

  % Filter
  for ink=1:28
   Lt(ink)=10*log10(sum((10.^(XdBn/10)).*(abs(H(ink,:).^2))));
  end 

  % Calculate Loudness
  [Nseg,Ns_seg,err] = DIN45631(Lt,'f');
  N(k) = Nseg;
  rmsvals(k) = rms(seg);
  %N2rms(k) = N(k)/rmsvals(k);
  for l=1:Nbrks
    N_spec_brks(l) = mean(Ns_seg(b1(l):b2(l)));
  end
  N_spec(:,k) = N_spec_brks';
end


%figure;
%subplot(311); plot([1:length(sp)]/Fs,sp);
%nz = length(find(N~=0));
%nz2 = length(find(N2rms~=0));
%subplot(312); plot((length(sp)/Fs)*[1:length(N)]/length(N),N);
%hold on; plot((length(sp)/Fs)*[1:length(N)]/length(N),N2rms,'r');
%legend(sprintf('%f',sum(N)/nz),sprintf('%f',sum(N2rms)/nz2))
%subplot(313); plot((length(sp)/Fs)*[1:length(N)]/length(N),rmsvals,'k');
%legend(sprintf('%f',mean(rmsvals)))

