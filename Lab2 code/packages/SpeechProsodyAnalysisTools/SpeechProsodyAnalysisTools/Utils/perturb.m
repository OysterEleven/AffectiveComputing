function [PF,PQ] = perturb(x,K)
%
% Finds the perturbation factor and perturbation quotient
% of a sequence
%
% Usage: [PF,PQ] = perturb(x,K);
%
%	x: Input signal.
%	K: Number of samples in the window used to find local mean in 
%          perturbation quotient. If K is a vector with Nk entries,
%          PQ returns Nk values for each value of K. Each K value must be 
%          an odd number.
%      	   
%      PF: Perturbation factor. 
%      PQ: Perturbation quotient.
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


if (~isempty(find(rem(K,2) == 0))) error('K must be an odd integer'); end
if (max(K) > length(x)) error('K must be less than or equal to the length of x'); end

%% Perturbation Factor
PF = mean(abs(diff(x,1) ./ x(2:end))); 

%% Perturbation Quotient for each value in K
for k=1:length(K)
  Kk = K(k);
  mu = filter(ones(1,Kk)/Kk,1,x);
  PQ(k) = mean(abs((x(1+(Kk-1)/2:end-(Kk-1)/2) - mu(Kk:end))./mu(Kk:end)));
end
