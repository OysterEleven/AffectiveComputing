function [intervals,disscurve] = timbre2disscurve(freq,amp,intervals)
%
% Calculates a dissonance curve as a function of interval for a given 
% spectral pattern of harmonics (specified by their location and amplitudes) 
% using the model described in 
%    - W. A. Sethares, "Tuning, Timbre, Spectrum, Scale." Springer, 2004
%
%
% [I,D] = timbre2disscurve(freq,amp,interv)
%
% Inputs:
%      freq:  Nx1 vector containing location of harmonics
%       amp:  Nx1 vector containing the amplitude of the harmonics in freq
% intervals:  (optional) vector argument specifying the interval values
%             at which the perceptual dissonance curve should be evaluated
%
% Outputs:
%         I:  Kx1 vector with the interval values at which the perceptual
%             dissonance curve is evaluated. If interv is given as an input
%             parameter, then I=interv.
%         D:  Kx1 vector with the dissonance values at the intervals listed
%             in I.
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
  intervals = [1:0.01:2];
end

freq = freq(:)';
amp = amp(:)';

disscurve=[]; 

% 
% call function dissmeasure for each interval 
%                        
%for alpha=1+inc:inc:range 
for a = 1:length(intervals)
  alpha = intervals(a);
  f=[freq alpha*freq]; 
  a=[amp, amp]; 
  d=dissmeasure(f, a); 
  disscurve=[disscurve d]; 
end 

%plot(1:inc:range,diss) 
