function specparm = psp(y);
% 
% This function returns the Parabolic Spectral Parameter (PSP)
% associated with one cycle of the derivative of the glottal source 
% waveform
%
% Usage:  prm = psp(x)
%
%	x:  input signal assumed to correspond to one cycle 
%	    of the glottal source (derivative) from the instant
%	    of glottal opening to glottal closure. See GLOTTMARKS
%           for details on how to extract this information from an
%           input speech signal.
%     prm:  the parabolic spectral parameter
%
% Refs: "Parabolic Spectral Parameter -- A New Method for Quantification
%        of the Glottal Flow." Alku, P., Strik, H. and Vilkman, E.
%        Speech Communication, 1997, (22), 67-79
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



y = y(:)';
y = [y y(1)];
y = y-min(y);
y2 = y/sqrt(sum(y.^2));

% Find spectrum of glottal (derivative) waveform
Y = fft(y2,1024);
Y = 20*log(abs(Y));
Y = Y(1:512);


% Find parabola parameters for glottal derivative waveform
N = 3;
NE = [];
[a,b,ne] = parabparam(Y,N);
NE = [NE ne];
MaxIter = 512/2;
while ((ne<0.01) & (N < MaxIter))
  N = N+1;
  [a,b,ne] = parabparam(Y,N);
  NE = [NE ne];
end
a_wvf = a;

%% Repeat for DC signal

% Create DC signal of period equal to the length of y
% (y is assumed to consist of only one period)
% and find its spectrum
dc = ones(1,length(y));
dc = dc/sqrt(sum(dc));
DC = fft(dc,1024);
zeroind = find(DC == 0);
DC(zeroind) = 0.0001;
DC = 20*log(abs(DC));
DC = DC(1:512);

% Find optimal parabolic parameter for DC signal
N = 3;
NE = [];
[a,b,ne] = parabparam(DC,N);
NE = [NE ne];
while (ne<0.015)
  N = N+1;
  [a,b,ne] = parabparam(DC,N);
  NE = [NE ne];
end
a_max = a;

% Normalize PSP of the waveform by PSP of DC signal
specparm = a_wvf/a_max;


function [a,b,NE] = parabparam(Y,N);

Y=Y(:)';
kN = [0:N-1].^2;
XN = Y(1:N);
a = (N*sum(XN.*kN) - sum(XN)*sum(kN))/(N*sum(kN.^2)-sum(kN)^2);
b = sum(XN-a*kN)/N;

NE = sum((XN-a*kN-b).^2)/sum(XN.^2);


