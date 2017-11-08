function [Timevals,Intervals,Disscurves,IntrinsicDiss] = dissgram(sp,Fs,MaxNumHarm,Intervals,F0,Tind) 
%
% [Timevals,Intervals,Disscurves,IntrinsicDiss] 
%                = dissgram(sp,Fs,MaxNumHarm,Intervals,F0,Tind) 
%
% Calculates a dissonance diagram as a function of time and frequency
% intervals. Each column of the diagram is the perceptual dissonance
% curve for the spectral pattern obtained for a short-term window frame
% of speech. See W. A. Sethares, "Tuning, Timbre, Spectrum, Scale." 
% Springer, 2004; and R. Fernandez "A Computational Model for the 
% Recognition of Affect from Speech" MIT Ph.D. Thesis 2003. 
%
% Inputs:
%	   sp:	speech
%          Fs:	sampling frequency
%  MaxNumHarm:  maximum number of harmonics to include in dissonance
%		calculation (Default: 6)
%   Intervals:  Nx1 array containing frequency interval values
%          F0:	pitch contour (if available from a call to getF0; 
%               to compute, default to [])_ 
%        Tind:  corresponding frame information for F0 (if available;
%               skip this argument to compute). See getF0 for more
%               help on the last two arguments.
%
%  Outputs:
%    Timevals:  Mx1 array containing time value of mth speech frame
%   Intervals:  Nx1 array containing frequency interval values 
%               (Same as input argument)
%  Disscurves:  NxM array where the nth column contains the perceptual
%               dissonance curve for the spectral pattern obtained for
%               the mth speech frame.
% IntrinsicDiss: Mx1 array with the intrinsic dissonance value for each
%                speech frame
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


if nargin < 3 | isempty(MaxNumHarm); MaxNumHarm = 6; end
if nargin < 4 | isempty(Intervals); Intervals = [1:0.01:2.5]; end
if nargin < 6 | isempty(F0) | isempty(Tind)
  [F0,str,Tind,wf] = getF0(sp,Fs,[],[],450);
end

vind = find(F0 ~= 0);
NumFrames = length(F0);
Nfft = 512;

IntrinsicDiss = zeros(1,NumFrames);
for k=1:length(vind)
  % Get voiced speech segment
  seg = sp(Tind(vind(k),1):Tind(vind(k),2));

  % Find PSD
  Pxx = psd(seg,Nfft,Fs);

  % Find location and amplitude of first 'MaxNumHarm' harmonics
  f0bin = round(F0(vind(k))*Nfft/Fs);
  offset = round(0.75*f0bin);
  spacing = round(0.55*f0bin);
  endset = offset + (MaxNumHarm+2)*spacing;
  [l,v] = mypickpeak(Pxx,MaxNumHarm,spacing,[],offset,endset);

  % Normalize so that highest harmonic is 1, and convert to Hz
  freqamps = v/max(v);
  freqloc = (l-1)*Fs/Nfft;     %Subtract 1 b/c DC is first entry


  % Find dissonance curve for this harmonic pattern
  [Intervals,disscrv] = timbre2disscurve(freqloc,freqamps,Intervals);
 
  % Testing stuff 
  %subplot(211); plot(Pxx)
  %subplot(212); plot(intervals,disscrv)
  %pause
  Disscurves(:,vind(k)) = disscrv';
  IntrinsicDiss(vind(k)) = dissmeasure(freqloc',freqamps');
end
[nr,nc] = size(Disscurves);
Disscurves = [Disscurves zeros(nr,NumFrames-vind(end))];


Timevals = mean(Tind')/Fs;

%figure
%subplot(311)
%imagesc(Timevals,Intervals,sqrt(sqrt(Disscurves)))
%colormap hot
%axis xy
%subplot(312)
%plot(Timevals,F0,'r*')
%subplot(313)
%plot([1:length(sp)]/Fs,sp,'r')

