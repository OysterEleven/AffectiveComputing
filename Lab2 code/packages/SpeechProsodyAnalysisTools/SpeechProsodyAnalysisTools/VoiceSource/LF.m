function [g,J] = LF(theta,x,Ne,N0);
%
% Calculates the glottal pulse derivative waveform (and Jacobian) of
% the parametric Liljencrants-Fant model.
% 
% Usage:
%     g = LF(theta,x,Ne,N0)  or     g = LF(theta,[],Ne,N0)
% [g,J] = LF(theta,x,Ne,N0)  or [g,J] = LF(theta,[],Ne,N0)
%
%   theta: vector containing the four parameters of the piecewise LF
%          model [alpha,Tp,Ee,beta], where alpha and Tp are parameters in
%          y=E0*exp(alpha*y)*sin(pi*y/Tp) for 0<y<Ne,
%          and Ee and beta are parameters for Ne<y<N0 given by the return phase
%          y=(Ee/(1-exp(-beta*(N0-Ne))))*(exp(-beta*(y-Ne))-exp(-beta*(N0-Ne))).
%          N0 is the length (period) of the signal, and E0 is a constant
%          determined by E0 = -Ee / (exp(alpha*Ne)*sin(pi*Ne/Tp));
%       x: vector containing time samples at which the model is evaluated
%          Default is 1:1:N0 (i.e., sampling period is assumed to be 1 w.l.o.g.)
%      Ne: instant of maximum excitation (in samples): 1<Ne<N0
%      N0: fundamental period (in samples)
%
%       g: N0x1 vector containing one period of the glottal pulse derivative
%       J: N0x4 vector with the i-th column containing the Jacobian of g with 
%          respect to the i-th parameter in alpha at each time sample. J is an
%	   optional argument (used by fitLF during the optimization to evaluate
%          gradients given a certain parameter estimate)
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




if isempty(x)
  x=[1:1:N0];
end

alpha = theta(1);
Tp = theta(2);
Ee = theta(3);
beta = theta(4);

%  Calculate the values of the LF model
x1 = x(1:1:Ne);
x2 = x(Ne+1:1:N0);
E0 = -Ee / (exp(alpha*Ne)*sin(pi*Ne/Tp));
g1 = E0*exp(alpha*x1).*sin(pi*x1/Tp);
g2 = -(Ee/(1-exp(-beta*(N0-Ne))))*...
      (exp(-beta*(x2-Ne))-exp(-beta*(N0-Ne)));
g=[g1 g2]';

% Calculate Jacobian if requested
if nargout > 1

  % w.r.t. alpha
  d_alpha1 = -(Ee/sin(pi*Ne/Tp)).*exp(alpha*(x1-Ne))...
             .*(x1-Ne).*sin(pi*x1/Tp);
  d_alpha2 = 0*x2;
  d_alpha = [d_alpha1 d_alpha2];

  % w.r.t. Tp
  d_Tp1 = (x1.*cos(pi*x1/Tp)*sin(pi*Ne/Tp) - Ne*sin(pi*x1/Tp)*cos(pi*Ne/Tp));
  d_Tp1 = (Ee*pi*d_Tp1.*exp(alpha*(x1-Ne)))/((Tp*sin(pi*Ne/Tp))^2);
  d_Tp2 = 0*x2;
  d_Tp = [d_Tp1 d_Tp2]; 
 
  % w.r.t. Ee
  d_Ee1 = -exp(alpha*(x1-Ne)).*sin(pi*x1/Tp)/sin(pi*Ne/Tp);
  d_Ee2 = (exp(-beta*(N0-Ne))-exp(-beta*(x2-Ne)))/(1-exp(-beta*(N0-Ne)));
  d_Ee = [d_Ee1 d_Ee2];

  % w.r.t. beta
  d_beta1 = 0*x1;
  d_beta2 = (x2-Ne).*exp(-beta*(x2-Ne))-x2.*exp(-beta*(x2-2*Ne+N0))...
          -(N0-Ne)*exp(-beta*(N0-Ne));
  d_beta2 = Ee*d_beta2/((1-exp(-beta*(N0-Ne)))^2);
  d_beta = [d_beta1 d_beta2];

  % Create Jacobian
  J = [d_alpha; d_Tp; d_Ee; d_beta]';
end

