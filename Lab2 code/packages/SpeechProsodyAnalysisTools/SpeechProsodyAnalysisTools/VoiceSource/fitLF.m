function [theta,gh,eop,ecp] = fitLF(yn,Ne,verbos)
%
% Fits the Liljencrants-Fant model to an estimate of the glottal
% pulse derivative (i.e., obtained via inverse-filtering,etc.)
% This algorithm uses the optimization toolbox to solve a non-linear
% minimization problem with bounded constraints on the parameters.
%
% Usage: theta = fitLF(y)         or     theta = fitLF(y,Ne)
%    [theta,g] = fitLF(y)         or [theta,g] = fitLF(y,Ne)
%    [theta,g,eop,ocp] = fitLF(y) or [theta,g,eop,ecp] = fitLF(y,Ne)
%
% 	y: one period of a glottal pulse derivative
%      Ne: the location of the instant of maximum excitation.
%          The default is the location of the minimum of y.
%          Note that this parameter is not optimized. For reliable
%          identification of Ne, see EXCITATION_INSTANTS
%  verbos: Level of verbosity: 'off', 'iter', 'notify', 'final' (default).
%
%   theta: vector containing the four parameters of the piecewise LF 
%          model [alpha,Tp,Ee,beta], where alpha and Tp are parameters in
%    	   y=E0*exp(alpha*y)*sin(pi*y/Tp) for 0<y<Ne, 
%          and Ee and beta are parameters for Ne<y<N0 given by the return phase
%          y=(Ee/(1-exp(-beta*(N0-Ne))))*(exp(-beta*(y-Ne))-exp(-beta*(N0-Ne))).
%          N0 is the length (period) of the signal, and E0 is a constant
%          determined by E0 = -Ee / (exp(alpha*Ne)*sin(pi*Ne/Tp));
%
%       g: optional output argument returning the synthetic LF 
%	   pulse derivative for the estimated parameters, with Ne and N0
%          as in y.
%
%     eop: optional output argument returning the mean squared error between
%          y and the LF fit throughout the open phase only;
%
%     ocp: optional output argument returning the mean squared error between
%          y and the LF fit throughout the closed phase only.
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



yn = yn(:);
N0 = length(yn);
if nargin < 2 | isempty(Ne)
  n1 = round(0.15*N0); n2 = round(0.75*N0);
  [miny,Ne] = min(yn(n1:n2));
  Ne = Ne + n1 - 1;
end


if Ne >= N0
  error('Ne parameter must be smaller than the length of the signal')
end

if nargin < 3 | isempty(strcmp(verbos,{'off', 'iter', 'notify', 'final'}))
  verbos = 'final';
end

% Initial parameter estimate
theta_0 = [0.0001 .75*Ne -min(yn) 3];

% Set options and bounds for optimization
options = optimset('Jacobian','on','LargeScale','on','Display',verbos,'TolX',1e-7);
lb = [0 0.51*Ne 0 0];
ub = [1 0.95*Ne Inf 5];
 
% Optimize
[theta,resnorm,resdl,exitflag,out] = lsqcurvefit(@LF,theta_0,[1:N0],yn,lb,ub,options,Ne,N0);

% Synthesize pulse corresponding to estimated parameters
if nargout > 1
  gh = LF(theta,[],Ne,N0);
end

if nargout > 2
  eop = mean((gh(1:Ne)-yn(1:Ne)).^2);
  ecp = mean((gh(Ne+1:N0)-yn(Ne+1:N0)).^2);
end
