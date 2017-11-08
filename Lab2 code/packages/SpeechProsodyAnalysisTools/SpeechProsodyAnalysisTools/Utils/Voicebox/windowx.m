function [v,bw] = windowx(nbw,n,b,z)
%WINDOWX Generate ordinates for WINDOW()

% return value is a vector equal to (A:B:C)'/D
% z = 1 if window goes to zero, else z=0 and modes 0 & 4 become 1 & 5
%
%						     Window Function
%
% b1 b2 b3	 A   B   C    D   	----====::::####<<<<$>>>>####::::====----
%
% 0  0  0	1-m  2	m-1  n+1  	-   *   *   *   *   *   *   *   *   *   -                 
% 0  0  1	1-m  2	m-1  n-1  	*   *   *   *   *   *   *   *   *   *   *                 
% 0  1  0	1-m  2	m-1   n	  	  *   *   *   *   *   *   *   *   *   *
% 0  1  1	 -m  2	m-2   n  	*   *   *   *   *   *   *   *   *   *                     
% 1  0  0	 0   1	m-1   n  	                    *   *   *   *   *   -
% 1  0  1	 0   1	m-1  n-1  	                    *   *   *   *   *   *
%
%	Mike Brookes, Imperial College April 1995


%      Copyright (C) Mike Brookes 2002
%
%      Last modified Mon Jan  7 13:59:37 2002
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

b1=b(1);b2=b(2);b3=b(3);

kb=2-b1;
if n<1,
  bw=n;
  n=kb*nbw/n+1-b2;
  kd=n+b2-1;
  m=floor(n);
else
  m=floor(n);
  if n>m, z=0; end;
  kd=n+(1-b2)*(z*(1-b3)*kb-1);
  bw=kb*nbw/kd;end;

v = ((b1-1)*(m-1)-b2*b3:kb:m-1-b2*b3)'/kd;

