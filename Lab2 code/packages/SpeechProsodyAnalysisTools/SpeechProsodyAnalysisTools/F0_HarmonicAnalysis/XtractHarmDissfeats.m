function harmdissfeats = XtractHarmDissfeats(dissgramrecord,N0,N1)
%
% harmdissfeats = XtractHarmDissfeats(dissgramrecord,N0,N1)
%
% dissgramrecord:  a structure containing the fields
%		    . time
%		    . intervals
%		    . dissgram
%		    . IntDiss
%		    . F0
%		    . F0frames
%          N0,N1:  endpoints (in terms of speech samples) of analysis window	
%


%
% Copyright (C) 2006  Raul Fernandez
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


%% Find analysis window endpoints
framecenters = mean(dissgramrecord.F0frames')';
if nargin < 2
 N0 = 1; N1 = framecenters(end);
end
t0 = find(framecenters >= N0);
t1 = find(framecenters <= N1);
if isempty(t1) | isempty(t0)
 disp('Analysis interval too small. Zeroing all features. Exiting')
 harmdissfeats = zeros(1,14);
 return;
else
 t0=t0(1); t1=t1(end);
end

 
%% Trim the corresponding variables related to the time axis
T = dissgramrecord.time;
T = T(t0:t1);
I = dissgramrecord.intervals;
D = dissgramrecord.dissgram;
D = D(:,t0:t1);
IntDiss = dissgramrecord.IntDiss;
IntDiss = IntDiss(t0:t1);

feat = local_disscurvefeats(T,I,D);
vind = find(sum(feat)~=0);

if ~isempty(vind)
  medfeat1 = median(feat(:,vind)');
  IntDissmed = median(IntDiss(vind));
  IntDissrange = max(IntDiss(vind)) - min(IntDiss(vind));
  harmdissfeats = [IntDissmed IntDissrange medfeat1];
else
  harmdissfeats = zeros(1,14);
end

%%%%% Local function   %%%%%%%%%

function feats = local_disscurvefeats(T,I,D)


vind = find(D(3,:)> 0);

Nf = 12;
feats = zeros(Nf,length(T));
for k=1:length(vind)
   %Select a disscurve slice (normalize so that max is 1?)
   dc = D(:,vind(k));

   % Area under curve
   dcarea = sum(dc)/length(dc);
   dcdiffarea = mean(abs(diff(dc)));

   % Find local valleys (min) and exclude unison
   % Retain two deepest valleys
   [l,v] = mypickpeak(-dc,15,2,-max(dc));
   if l(1) == 1; l(1)=[]; v(1)=[]; end
   vlminsort = sortrows([-v l]);
   if ~isempty(vlminsort)
     Valley1val = vlminsort(1,1);
     Valley1loc = I(vlminsort(1,2));

     Valley2val = Valley1val; Valley2loc = Valley1loc;
     if (size(vlminsort,1) > 1)
       Valley2val = vlminsort(2,1);
       Valley2loc = I(vlminsort(2,2));
     end

     Valleymuval = mean(vlminsort(:,1));
     Valleymuloc = mean(I(vlminsort(:,2)));
   else
     Valley1val = 0; Valley1loc = 0;
     Valley2val = Valley1val; Valle2loc = Valley1loc;
     Valleymuval = Valley1val; Valleymuloc = Valley1loc;
   end

   % Find local peaks
   % Retain two highest peaks
   [l,v] = mypickpeak(dc,15,2);
   if l(end) == length(I); l(end) = []; v(end)=[]; end
   vlmaxsort = sortrows([v l]);
   vlmaxsort(end,:) = [];
   if ~isempty(vlmaxsort)
     Peak1val = vlmaxsort(end,1);
     Peak1loc = I(vlmaxsort(end,2));

     Peak2val = Peak1val; Peak2loc = Peak1loc;
     if size(vlmaxsort,1) > 1
       Peak2val = vlmaxsort(end-1,1);
       Peak2loc = I(vlmaxsort(end-1,2));
     end

     Peakmuval = mean(vlmaxsort(:,1));
     Peakmuloc = mean(vlmaxsort(:,2));

   else
     [vm,lm] = max(dc);
     Peak1val = vm; Peak2val = vm; Peakmuval = vm;
     Peak1loc = I(lm); Peak2loc = Peak1loc; Peakmuloc = Peak1loc;
   end

   %dccurvefeats = [dcarea dcdiffarea minimavalmu minimalocmu maximavalmu maximalocmu ...
   %                maxminval maxminloc minminval minminloc maxmaxval maxmaxloc minmaxval minmaxloc];

   dccurvefeats = [dcarea dcdiffarea ...
                   Valley1val Valley1loc Valley2val Valley2loc Valleymuval ...
                   Peak1val Peak1loc Peak2val Peak2loc Peakmuval];


   feats(:,vind(k)) = dccurvefeats';
end

