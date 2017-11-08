function [loc,val] = mypickpeak(X,max_num_pks,min_pk_spc,min_thresh,n0,nT)
%
% [loc,val] = mypickpeak(X,max_num_pks,min_pk_spc,min_thresh,n0,nT)
%
%             X:   signal
%   max_num_pks:   maximum number of peaks (Default is 10) (if fewer are found, NaN values are omitted)
%    min_pk_spc:   minimum spacing (in samples) between adjacent peaks (Default is 2)
%    min_thresh:   minimum height for a peak to be valid  (Default is 0)
%            n0:   initial sample to search peaks over (Default is 1)
%            nT:   final sample to search peaks over (Default is the length of X)
%
%           loc:   max_num_pks x 1 vector (or smaller) containing the peak indices in
%                  ascending order
%           val:   vector of the same length as loc, containing the respective peak heights
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


X = X(:)';
nc = length(X);
if ((nargin < 2) | isempty(max_num_pks)) max_num_pks = 10; end
if ((nargin < 3) | isempty(min_pk_spc)) min_pk_spc = 2; end
if ((nargin < 4) | isempty(min_thresh)) min_thresh = 0; end
if ((nargin < 5) | isempty(n0)) n0=1; end
if ((nargin < 6) | isempty(nT)) nT=nc; end

if (n0<1); n0=1; end
if (nT>nc); nT=nc; end


% Find peaks. 
[pk_loc,pk_hgt] = findpeak(X(n0:nT),max_num_pks,min_pk_spc);
pk_loc = pk_loc+n0-1;

% Discard 'unfound' peaks
disc_ind = find(isnan(pk_hgt));
if ~isempty(disc_ind)
  pk_loc(disc_ind) = []; pk_hgt(disc_ind)=[];
end


% Get rid of peaks smaller than threshold
tiny_pks_ind = find(pk_hgt <= min_thresh);
if ~isempty(tiny_pks_ind)
  pk_loc(tiny_pks_ind) = []; pk_hgt(tiny_pks_ind)=[];
end

P_LH = sortrows([pk_loc pk_hgt],1);
loc = P_LH(:,1); val = P_LH(:,2);


