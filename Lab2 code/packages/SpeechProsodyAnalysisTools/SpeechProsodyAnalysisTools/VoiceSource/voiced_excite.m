function vexcit_struct = voiced_excite(x,Fs)
%
% This function returns a voiced excitation structure containing a 
% segmentation of a speech waveform into non-overlapping voiced segments 
% and information about the instants of voiced excitation for each segment. 
% (By 'segment' it is meant a portion or fragment of speech, not a 
% phonological or even phonetic segment.)
%
% Usage:
%          vexcit = voiced_excite(x,Fs)
%
% Inputs:
%	x: input signal (Nx samples in length)
%      Fs: sampling frequency
%
% Output:
%  vexcit: a structure containing the following fields:
%            -num_segs:    the number of voiced segments found in x
%            -excitepulse: a pulse train corresponding to instants of voiced 
%                          excitations (an Nx length row vector)
%            -exciteinds:  a 1xNvs cell array, where the kth cell contains 
%                          indices of the voiced excitations within the kth 
%                          segment
%            -f0segs:      a 1xNvs cell array, where the kth cell contains
%                          the smoothed-out f0 contour corresponding to the 
%                          k- th segment
%
% Note: Because the information returned by voiced_excite is fairly 
% computationally intensive, and its output can be used by several other
% functions in this toolbox, it is convenient to obtain a voiced excitation 
% structure through this function first before invoking other functions. 
% To save cycles (and time), pass it as an optional parameter when 
% subsequently calling other functions that would internally recompute these 
% values if not specified. See, e.g., JITT_SHIMM.
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



Nx = length(x);


%%% Get pitch contour, strenghts, etc.
[F0,strength,T_ind] = getF0(x,Fs);


%%% Find voiced frames
voiced_ind = find(F0~=0);
voicing_breaks=find(diff(voiced_ind)>1);
VI(1,1) = voiced_ind(1);
k=0;
if ~isempty(voicing_breaks)
  num_voiced_segm = length(voicing_breaks);
  for k=1:num_voiced_segm
     VI(k,2) = voiced_ind(voicing_breaks(k));
     VI(k+1,1) = voiced_ind(voicing_breaks(k)+1);
  end
end
VI(k+1,2)=voiced_ind(end);


%%% Get rid of "spurious" voiced regions (i.e., regions shorter than N frames)
N = 3;
framediff = VI(:,2)-VI(:,1);
shortsegs = find(framediff < N);
if ~isempty(shortsegs) VI(shortsegs,:) = []; end


%%% Find the excitation instant for each voiced segment
Nvs = size(VI,1);
excitations = zeros(1,Nx);
p2pinterval=[];
time_inds = [];
for k=1:Nvs
  voiced_seg = x(T_ind(VI(k,1)+1,1):T_ind(VI(k,2)-1,2));
  f0segs{1,k} = F0(VI(k,1):VI(k,2));
  meanF0k = mean(F0(VI(k,1):VI(k,2)));
  [pks,slps,res,sig]=excitation_instants(voiced_seg,Fs,meanF0k);
  excitations(pks+T_ind(VI(k,1)+1,1)) = ones(1,length(pks));
  excitations_inds{1,k} = pks+T_ind(VI(k,1)+1,1);
end 

%% Bundle relevant variables into the voiced excitation structure
vexcit_struct.num_segs=Nvs;
vexcit_struct.excitepulse = excitations;
vexcit_struct.exciteinds = excitations_inds;
vexcit_struct.f0segs = f0segs;


