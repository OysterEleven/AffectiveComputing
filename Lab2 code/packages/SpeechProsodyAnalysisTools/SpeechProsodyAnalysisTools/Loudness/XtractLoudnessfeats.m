function loudnessfeats = XtractLoudnessfeats(in_Nloud,in_RMSCurve,in_Tind,in_Barks,in_N_spec,N0,N1)
%
% loudnessfeats = XtractLoudnessfeats(Nloud,RMSCurve,Tind,Barks,N_spec,N0,N1)
%
% Extracts a set of features from a loudgram.
%
% Inputs:
%           All but the last two inputs to this function are the outputs 
%           of loudnessgram (see LOUDNESSGRAM for  explanation). N0 and N1 
%           indicate the first and last speech samples defining the segment
%           over which the features are extracted. If either is empty, the 
%           analysis is done over the  entire contour.
%
% Output:
% loudnessfeats is a 20x1 vector containing statistics derived from the
% loudness contour; namely:
%    (1)   : the area under the loudness curve
%    (2-4) : the 25-, 50- and 75-percentile of the absolute loudness curve
%    (5-7) : the 25-, 50- and 75-percentile of the RMS curve
%    (8-20): the mean value of the specific loudness for Barks 5 through 17 
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



% Get indices of analysis region
if (nargin < 6 | isempty(N0) | isempty(N1))
  N0 = 1;
  N1 = length(in_Tind);
else %Convert speech samples to analysis frames in the loudgrams
  framecenters = mean(in_Tind')';
  T0 = find(framecenters>=N0);
  if isempty(T0)
    error('N0 lies beyond end of speech segment');
  end
  N0 = T0(1);
  T1 = find(framecenters<=N1);
  if (isempty(T1))
    error('No analysis data points found before N1');
  end
  N1 = T1(end);
end 
  	
% Extract loudness measures over the analysis interval 
Nloud    = in_Nloud(N0:N1);
RMScurve = in_RMSCurve(N0:N1);
Barks    = in_Barks;
N_spec   = in_N_spec(:,N0:N1);



nzind = find(Nloud ~= 0);
Nnz = Nloud(nzind);
if isempty(nzind)
  loudnessfeats = zeros(1,12);
  return;
end

% Area
loudnessarea = mean(Nnz);


% Percentiles of absolute loudness curve
p = [25 50 75];
for l=1:length(p)
  pth = prctile(Nnz,p(l));
  Nprctilemean(l) = mean(Nnz(find(Nnz > pth)));
end

% Percentiles of rms curve
RMSnz = RMScurve(nzind);
for l=1:length(p)
  pth = prctile(RMSnz,p(l));
  RMSprctilemean(l) = mean(RMSnz(find(RMSnz > pth)));
end
 
% Get median of the specific loudness for barks 5 to 17 
%(freq range 510Hz-3.7kHz)
barkrange = [4.25:1:17.25];
for l=1:length(barkrange)-1
  bark1ind = find(Barks==barkrange(l));
  bark2ind = find(Barks==barkrange(l+1));
  barkenergy(l)=mean(mean(N_spec(bark1ind:bark2ind,nzind),1));
end 


% Loudness feats
loudnessfeats = [loudnessarea Nprctilemean RMSprctilemean barkenergy];



