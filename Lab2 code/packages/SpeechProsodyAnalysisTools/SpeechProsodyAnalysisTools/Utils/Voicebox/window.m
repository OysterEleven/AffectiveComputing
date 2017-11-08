function [w,bw0] = window(wtype,n,mode,p)
%WINDOW Generate a standard windowing function
%	[W,BW] = WINDOW(TYPE,N or BW,MODE,P) calculates a window function
%
%	TYPE	is one of:
%			'blackman'
%			'kaiser'	with parameter P (often called beta) (default P=8)
%			'gaussian'	truncated at P std deviations (default P=3)
%			'hamming'
%			'hanning'
%			'harris3'	3-term balckman-harris with 67dB sidelobes
%			'harris4'	4-term balckman-harris with 92dB sidelobes
%			'rectangle'
%			'triangle'
%
%	N or BW	is either the number of points to generate or else the desired
%		bandwidth divided by the sample frequency. The bandwidth used is
%		the 6dB bandwidth which gives the required separation of unit sinewaves
%		for the Fourier transform to produce a dip between them. If N>1 is
%		a non-integer, then FLOOR(N) points will be generated.
%
%	MODE	is used to determine the exact placing of the generated points within
%		the selected window function. The six options available are illustrated
%		below. Mode 0 is the default.
%						     Window Function
%
%		Mode Period     ----====::::####<<<<$>>>>####::::====----
%
%		 0	N+1,N-1		-   *   *   *   *   *   *   *   *   *   -                 
%		 1	  N-1		*   *   *   *   *   *   *   *   *   *   *                 
%		 2	   N		  *   *   *   *   *   *   *   *   *   *                     
%		 3	   N		*   *   *   *   *   *   *   *   *   *
%		 4	2N,2N-2		                    *   *   *   *   *   -
%		 5	 2N-2		                    *   *   *   *   *   *
%
%	Notes:	a) Modes 0 and 4 omit the end point of the window if it is zero. For
%		   window functions that don't go to zero, these modes are the same as 1 and 5.
%
%		b) The central point of the window is included in modes 0, 1 and 2 iff n is odd,
%		   is included in mode 3 iff n is even and is always included in modes 4 and 5
%
%		c) If N is a non-integer, then it will still determine the window function but
%		   only FLOOR(N) points will be generated.
%
%		d) All windows are scaled to have a peak value of 1.
%
%	P	parameter needed for some windows
%
%	BW	returns the 6dB bandwidth of the window. Multiply
%		by the sample frequency to get the bandwidth in Hz.
%		Bandwidth values are only accurate to a couple of significant figures.

%	Mike Brookes, Imperial College April 1995
% Ref: Proc IEEE Vol 66 No 1 Jan 1978 pp51-83

% Test1: f=20*log10(cliplow(abs(fft([window('xxx',100,1); zeros(700,1)])),0.0001)); plot(f(1:100)-f(1))
% Test2: for i=0:5,[w,bw]=window('gaussian',6,i);v(i+1,1:7)=[bw;w]';end;v


%      Copyright (C) Mike Brookes 2002
%
%      Last modified Mon Jan  7 13:59:27 2002
%
%   VOICEBOX is a MATLAB toolbox for speech processing. Home page is at
%   http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3, mode=0; end;
b = rem(floor(mode*pow2(-2:0)),2);
wtype = lower(wtype);

if strcmp(wtype,'hanning'),
  [v,bw]=windowx(1,n,b,1);
  w = 0.5+0.5*cos(pi*v);

elseif strcmp(wtype,'rectangle'),
  [v,bw]=windowx(0.605,n,b,0);
  w = ones(size(v));

elseif strcmp(wtype,'triangle'),
  [v,bw]=windowx(0.89,n,b,1);
  w = 1-abs(v);

elseif strcmp(wtype,'gaussian'),
  if nargin<4, p=3; end;
  [v,bw]=windowx(0.36*p,n,b,0);
  w=exp(-0.5*p^2*(v.*v));

elseif strcmp(wtype,'kaiser'),
  if nargin<4, p=8; end;
  [v,bw]=windowx(0.12*p,n,b,0);
  w=besseli(0,p*sqrt(1-v.^2))/besseli(0,p);

elseif strcmp(wtype,'hamming'),
  [v,bw]=windowx(0.905,n,b,0);
  w = 0.54+0.46*cos(pi*v);

elseif strcmp(wtype,'blackman'),
  [v,bw]=windowx(1.175,n,b,1);
  w = 0.42+0.5*cos(pi*v) + 0.08*cos(2*pi*v);

elseif strcmp(wtype,'harris3'),
  [v,bw]=windowx(0.905,n,b,0);
  w = 0.42323 + 0.49755*cos(pi*v) + 0.07922*cos(2*pi*v);

elseif strcmp(wtype,'harris4'),
  [v,bw]=windowx(1.36,n,b,0);
  w = 0.35875 + 0.48829*cos(pi*v) + 0.14128*cos(2*pi*v) + 0.01168*cos(3*pi*v);

end;
if nargout>1, bw0=bw;end;


