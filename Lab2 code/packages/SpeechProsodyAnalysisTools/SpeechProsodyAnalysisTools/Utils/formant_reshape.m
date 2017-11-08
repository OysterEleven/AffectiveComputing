function ar_out = formant_reshape(ar_in,Fs,lo_thresh,hi_thresh)
%
% ar_out = formant_reshape(ar_in,Fs,lo_thresh,hi_thresh)
%
% This function "fixes" an AR polynomial by discarding the roots that
% correspond to 
%   (i)  formants with a center frequency below lo_thresh (default: 50Hz)
%        and above hi_thresh (default:0.95*Fnyquist)
%   (ii) poles on the positive real axis
% The gain is normalized so that the square norm of the coefficients in 
% the original AR polynomial is preserved.
% 
% Usage: ar_out = formant_reshape(ar_in,Fs,lo_thresh,hi_thresh)
%
%         ar_in: original AR polynomial
%            Fs: sampling frequency
%     lo_thresh: center frequency lower threshold
%     hi_thresh: center frequency upper threshold
%        ar_out: modified AR polynomial
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


if (nargin < 2) error('Insufficient number of arguments'); end
if ((nargin < 3) | isempty(lo_thresh)) lo_thresh=50; end
if ((nargin < 4) | isempty(hi_thresh)) hi_thresh=0.95*Fs/2; end

% Get roots of polynomial
arzz = roots(ar_in);

% Get (normalized) frequencies and bandwidths from roots
f=angle(arzz)*0.5/pi;
b=-log(abs(arzz))/pi;

% Find negative real roots and save for later
pos_roots_ind = find(f == 0.5);
pos_roots=[];
if ~isempty(pos_roots_ind)
  f(pos_roots_ind) = [];
  b(pos_roots_ind) = [];
  pos_roots = arzz(pos_roots_ind);
  arzz(pos_roots_ind) = [];
end


% Discard positive real roots (0Hz)
zero_roots = find(f == 0);
if ~isempty(zero_roots)
  f(zero_roots) = [];
  b(zero_roots) = [];
  arzz(zero_roots) = [];
end

% Discard roots with center freqs. below freq_thresh
freq_ind = find((abs(f*Fs) < lo_thresh) | (abs(f*Fs) > hi_thresh));
if ~isempty(freq_ind)
  f(freq_ind) = []; 
  b(freq_ind) = []; 
  arzz(freq_ind)=[]; 
end


% Discard roots with bandwidths above bwidth_thresh
%bw_ind = find(b*Fs > bwidth_thresh);
%if ~isempty(bw_ind)
%  f(bw_ind) = []; 
%  b(bw_ind) = []; 
%  arzz(freq_ind) = [];
%end

arzz=[arzz; pos_roots];

% Resulting AR polynomial
ar2 = poly(arzz);

% Normalize gains (equate energy in ar coefficients, after discarding
% roots, with energy in original ar coeffs.)
ar_out = ar2*(norm(ar_in)/norm(ar2));


