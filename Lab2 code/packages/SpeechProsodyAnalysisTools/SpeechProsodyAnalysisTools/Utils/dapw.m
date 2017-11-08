function a = dapw(x,p,N,MaxIter,K,W)
%
% This function finds a (frequency weighted) discrete all-pole model
% to a signal x. Unlike LP methods, the weighted discrete all-pole model 
% (DAPW) produces a better spectral fit for signals whose spectra can
% be best described by a relatively small set of discrete values (such
% as voiced speech). See "Discrete All-Pole Modeling" by A. El-Jaroudi
% and J. Makhoul (IEEE Trans. Signal Proc. 39(2), 1991).
%
% Usage: A = DAPW(x,P,N,MaxIter,Freqs,FreqWts) 
%
%        x = input signal
%	 P = model order (default is 12)
%        N = max. number of spectral peaks to fit (default is 40)
%  MaxIter = maximum number of iterations to run if algorithm
%	     doesn't converge earlier (default is 100)
%    Freqs = a 1xK vector containing the boundaries between adjacent
%	     frequency bands in increasing order. The last value corresponds
%            to the Nyquist frequency, and all other values are expressed as
%            a ratio thereof. For example, Freqs=[1/3 1/2 1]
%  FreqWts = The weights corresponding to each of the bands defined by Freqs.
%	     For instance, FreqWts=[1 2 1] weighs the error in the middle band by
%            twice the amount as in the adjacent bands.
%            If Freqs and FreqWts are omitted, all frequencies are weighted equally.
%
%        A = 1x(P+1) vector containing the coefficients of an all-pole model. 
%            The gain is factored into the coefficients, so A(1) need not equal 1.
%
% To skip arguments and assign them default values, use [] 
%      (e.g.,  A = DAPW(x,[],[],[],Freqs,FreqWts))
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
 

if ((nargin < 2) | isempty(p)) p = 12; end
if ((nargin < 3) | isempty(N)) N = 40; end
if ((nargin < 4) | isempty(MaxIter)) MaxIter = 100; end
if ((nargin < 5) | isempty(K))
   K = 1; W = 1;
end
if (nargin == 5) 
  error('You must enter Freqs and FreqWts pairs. Type "help dapw"')
end
if (min(diff(K)) < 0)
  error('Frequency band boundaries must be in increasing order. Type "help dapw"')
end
if (K(end) ~=1) 
  error('Last frequency boundary must correspond to Nyquist frequency. Type "help dapw"')
end
if (length(K) ~= length(W))
  error('Vectors Freqs and FreqWts must be of the same length')
end

x = x(:)';
M = length(x);

% Find signal spectrum
fftsize=max(512,2^nextpow2(M)); 
X = fft(x,fftsize);
P = abs(X).^2;

% Pick (at most) N peaks from the spectrum
[wn,Pw]=mypickpeak(P(1:fftsize/2+1),N,1);
zeroind=find(wn==0);
if ~isempty(zeroind); wn(zeroind)=[]; Pw(zeroind)=[]; end
N = length(wn);
w_Pw = sortrows([wn Pw],1);
wn = w_Pw(:,1);
Pw = w_Pw(:,2);


% Frequency weighting
wn_ratio = (wn-1)/(0.5*fftsize);
k1 = find(wn_ratio<=K(1));
Wgts(k1) = W(1);
for k=2:length(K)
   kk = find((wn_ratio > K(k-1)) & (wn_ratio<=K(k)));
   Wgts(kk) = W(k);
end
Wgts = Wgts'/sum(Wgts);

%Pw = Pw .* Wgts;

% Sample the spectrum (i.e., retain only the peaks) and i
% invert it to get the sampled autocorrelation fxn. R
Pwsampled = zeros(1,fftsize);
Pwsampled(wn) = Pw;
Pwsampled(fftsize/2+2:end)=fliplr(Pwsampled(2:fftsize/2));
R = real(ifft(Pwsampled));
Rminv = inv(toeplitz(R(1:p+1)));


% Initialize coefficients by doing standard LPC 
% (using R obtained above)
a = levinson(R,p);
ainit = a;
alpha = 1;
E = 1e+150;
%%%% Iterations begin here %%%%
k=1;
minalpha = 1e-20;
epsilon = 1e-10;
convergence = 0;

while (~convergence & (k <= MaxIter))
   % Evaluate frequency response of coeffs. a at freqs. wn
   A = freqz(a,1,pi*(wn-1)/(0.5*fftsize));
   Ainv=zeros(1,fftsize);
   Ainv(wn) = Wgts./A;
   Ainv(fftsize/2+2:end)=fliplr(conj(Ainv(2:fftsize/2)));

   h = real(fft(Ainv));
   h = h/h(1);
   anew = (1-alpha)*a +alpha*(Rminv*h(1:p+1)')';
   
   % Evaluate error
   Phat = abs(freqz(1,anew,pi*(wn-1)/(0.5*fftsize))).^2;
   E(k+1) = (sum(Pw./Phat - log(Pw./Phat))-1)/N;
   delta = E(end)-E(end-1);

   % If error increased in this iteration, do not leave the iteration
   % just yet; decrease alpha to force the error to go down
   while (delta > 0)
     alpha = 0.5*alpha;
     if alpha < minalpha;  %If alpha needs to be made smaller than minalpha
        convergence = 1    %to decrease the error, claim convergence since
        E(k)=[];           %the change in the error would be minimal anyways
        return;
     end
     anew = (1-alpha)*a +alpha*(Rminv*h(1:p+1)')';
     % Evaluate error
     Phat = abs(freqz(1,anew,pi*(wn-1)/(0.5*fftsize))).^2;
     E(k+1) = (sum(Pw./Phat - log(Pw./Phat))-1)/N;
     delta = E(end)-E(end-1);
   end
   
   % If relative change of delta is smaller than epsilon, claim convergence
   if (abs(delta)/E(end) <= epsilon) convergence=1; end
   a = anew;
   k = k+1;
end
Phat = abs(freqz(1,a,pi*(wn-1)/(0.5*fftsize))).^2;
G = sqrt(sum(Pw ./ Phat)/N);
a = a/G;

