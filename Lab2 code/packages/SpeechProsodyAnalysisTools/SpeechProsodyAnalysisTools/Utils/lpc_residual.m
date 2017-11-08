function [ar,e,k,res]=lpc_residual(signal,order,framedetails,windowtype,dc,fade)
%%%%%%%%%%%%%%%%%%% lpc_residual.m %%%%%%%%%%%%%%%%%%%
%
% Calculates the residual resulting from linear 
% predictive analysis on a speech signal. Performs
% the LPC, inverse-filtering the result from the 
% original signal and subtracting. 
%
%                 ------ usage ------
%
% [ar,e,k,res]=lpc_residual(signal,+order+,+framedetails+,
%			    +windowtype+,+dc+,+fade+); 
%
% where + denotes optional arguments.
%
%                 ------ inputs ------
%
% signal: obviously
% order: order of LPC analysis to be calculated (default: 10)
% framedetails: [framespacing framewidth skip]
%    Perform LPC analysis using frames of width 
%    framewidth, spaced framespacing samples apart,
%    skipping skip samples at the beginning.   
%    (default: [80 200 0])
% windowtype: OPTIONAL window type to apply in the LPC analysis. 
%    Type 'help window' for a list of acceptable window types and
%    their associated character strings. (default: 'hanning') 
% dc(nf,1): OPTIONAL column vector with as many rows as 
%    frames in the LPC analysis, each denoting the
%    dc component to be subtracted from the frame before 
%    inverse filtering (default: 0)
% fade: OPTIONAL. Autoregressive coefficients will be 
%    linearly interpolated for fade samples either side 
%    of frame boundaries (default: 0)
%
%                 ------ outputs ------
%
% ar(nf,order+1): the LPC autoregression coefficients, with ar(1)=1
% e(nf): the energy in the residual; sqrt(e) is 'filter gain'
% k(nf,2): first and last indices of each analysis interval
% res(nf): residual resulting from the LPC analysis
%
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



%%% input handling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if (nargin==1)
     order=10;
     framedetails=[80 200 0];
     windowtype='hanning';
     dc=0;
     fade=0;
elseif (nargin==2)
     framedetails=[80 200 0];
     windowtype='hanning';
     dc=0;
     fade=0;
elseif (nargin==3)
     windowtype='hanning';
     dc=0;
     fade=0;
elseif (nargin==4)
     dc=0;
     fade=0;
elseif (nargin==5)
     fade=0;
end%if
%
%
if ( (size(order)~=[1 1]) | (size(framedetails)~=[1 3]) | (isnumeric(windowtype)==1) | (size(dc)~=[1 1]) | (size(fade)~=[1 1]) )
     error('lpc_residual: Argument size error. Doublecheck the parameter list.')
end%if
%
%
%%% processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calculate LPC
%

[ar,e,k]=lpcauto(signal,order,framedetails);
%
% inverse filter 
%
res=lpcifilt(signal,ar,k(:,1),dc,fade);
res=res(:)';




