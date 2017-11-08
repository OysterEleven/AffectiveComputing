function [F1,T] = formant1a(x,Nw,Nskip,P,Fs,F1_th)
%
% Calculates the first formant from a speech segment
% using a sliding covariance LPC analysis method.
% 
%  [F1,T] = formant1(x,Nw,Nskip,P,Fs,F1_th)
%
% Inputs
%	x:	speech segment
%      Nw:	length of analysis (rectangular) window (in samples) 
%   Nskip:      number of samples to skip between adjacent frames
%       P:	order used in the LPC analysis
%      Fs:      sampling frequency
%   F1_th:	upper threshold for F1 [default: 1050]
%
% Outputs
%      F1:	track of first formant estimates
%       T:	the nth row contains the indices corresponding to the 
%		initial and final speech samples used to estimate the
%		first formant in the nth frame
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



if (nargin < 6)
  F1_th = 1050;
end

Nx = length(x);
x = x(:)';
%x = [zeros(1,P+1) x];
Nx = length(x);

% Perform short-term LPC analysis using the covariance method
nf = floor((Nx-Nw+Nskip)/Nskip);
T = [P+1:Nskip:nf*Nskip]';
T = [T T+Nw-1];
[ar,e,dc] = lpccovar(x,P,T);

% Get formants that fall below F1_th Hz.
[N,F,A,B] = lpcar2fm(ar,5);
Ffs = Fs*F;
NN = length(N);
F1=Ffs(:,1); 

outlierind=find((F1 > F1_th));
goodind = setdiff([1:NN],outlierind);
F1mf = medfilt1(F1,4);
if (~isempty(outlierind) & ~isempty(goodind))
  for k=1:length(outlierind)
    badind=outlierind(k);
    [goodneighbor,kn]=min(abs(goodind-badind));
    F1(badind) = F1mf(kn(1));
  end
end
if ~isempty(goodind); F1=F1mf; end
 

