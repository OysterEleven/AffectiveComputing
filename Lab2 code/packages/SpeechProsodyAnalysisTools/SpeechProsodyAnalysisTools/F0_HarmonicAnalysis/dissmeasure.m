function d=dissmeasure(f,amp) 
% 
% Calculates the intrinsic dissonance value for a spectral pattern of 
% harmonics specified by their location and amplitudes
% Ref: W. A. Sethares, "Tuning, Timbre, Spectrum, Scale." Springer, 2004 
%
% D = dissmeasure(F,A)
%
% Inputs:
%	F: location of harmonics
%	A: amplitude of the harmonics
%
% Output:
%	D: Intrinsic dissonance value of spectral pattern.
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

 
            
Dstar=0.24; S1=0.0207; S2=18.96; C1=5; 
C2=-5; A1=-3.51; A2=-5.75; 
N=length(f); D=0; 
[f,ind]=sort(f); 
ams=amp(ind); 
            
for i=2:N; 
            
  Fmin=f(1:N-i+1); 
  S=Dstar./(S1*Fmin+S2); 
  Fdif=f(i:N)-f(1:N-i+1); 
  a=ams(i:N).*ams(1:N-i+1); 
  Dnew=a.*(C1*exp(A1*S.*Fdif)+C2*exp(A2*S.*Fdif));
  D=D+Dnew*ones(size(Dnew))'; 
end 
            
d=D; 

