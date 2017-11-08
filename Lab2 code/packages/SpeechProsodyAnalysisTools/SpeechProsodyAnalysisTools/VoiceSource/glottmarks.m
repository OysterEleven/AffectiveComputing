function glottal_landmarks = glottmarks(x,Fs)
%
% Segments a speech waveform into phonatory cycles (in voiced regions) and 
% identifies "glottal landmarks" associated with the different phases of 
% each phonatory cycle.
%
% Usage:
%       glottal_inds = glottmarks(x,Fs)
%
% Inputs:
%      	       x: speech signal
%             Fs: sampling frequency
%
% Output:
%   glottal_inds: an Nx4 matrix, where N is the number of cycles found. Each 
%                 row consists of a series of indices for each phonatory 
%                 cycle, describing the following:
%                   - the instant of glottal opening
%                   - the instant of maximum excitation
%                   - the beginning of the closed phase
%                   - the end of the closed phase.
% 		  When two cycles are adjacent, the instant of glottal 
%                 opening is assumed to begin one sample after the end of 
%                 the closed phase of the previous cycle.
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




% Obtain excitation instants
[F0,strength,T_ind] = getF0(x,Fs,0.01);

% Separate into voiced segments; VI contains indices of voiced frames
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

% Discard segments containing fewer than N voiced frames
N = 2;
tinysegs = find ((VI(:,2) - VI(:,1)) < N-1);
if ~isempty(tinysegs); VI(tinysegs,:) = []; end
% Find the excitation instant for voiced segments only
% and build matrix with glottal landmarks: glottal opening,
% max. excitation, beginning of closed phase, and end of
% closed phase (when two cycles are adjacent, glottal opening 
% of a cycle is assumed one sample after end of closed phase 
% of previous cycle).  
Nvs = size(VI,1);
cp_intervals = [];
glottal_landmarks = [];
Nvs
for k=1:Nvs
     k
     index0 = T_ind(VI(k,1)+1,1); 
     indexT = T_ind(VI(k,2)-1,2);
     voiced_seg = x(index0:indexT);
     meanF0k = mean(F0(VI(k,1):VI(k,2)));
     [pks,slope,res,sig]=excitation_instants(voiced_seg,Fs,meanF0k);
     if length(pks > 1)
       cp = closed_phases_adaptive(voiced_seg,pks,12,Fs);
     else
       cp = [];
     end
     if ~isempty(cp) % Check to take care of degenerate case (only one true spike found)
       mark_pts = [cp(1:end-1,3)+ones(size(cp,1)-1,1) cp(2:end,:)];
       glottal_landmarks = [glottal_landmarks; mark_pts+index0*ones(size(mark_pts))];
     end
end

