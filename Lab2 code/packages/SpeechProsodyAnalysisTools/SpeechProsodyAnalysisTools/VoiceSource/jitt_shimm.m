function [jitter,shimmer,flag] = jitt_shimm(x,Fs,Kj,Ks,vexcit,N0,N1)
% 
% Usage [jitter,shimmer,flag] = jitt_shimm(x,Fs,Kj,Ks,vexcite,N0,N1)
%
% Returns jitter and shimmer measures quantified with the perturbation 
% factor (PF) and perturbation quotient (PQ) for the voiced segments of 
% a signal x.
%
% Inputs:
% 	x: input signal
%      Fs: sampling frequency
%      Kj: vector containing window lengths for evaluating the jitter PQ;
%          it must only contain odd integers. Default is 3. Use [] to bypass.
%      Ks: vector containing window lengths for evaluating the shimmer PQ;
%          it must only contain odd integers. Default is 3.
% vexcite: an optional argument containing a structure with 
%          information about the voiced excitations in x. 
%          See VOICED_EXCITE for details. Pass this argument if
%          already available from other processing to avoid recomputing.
%      N0: beginning of analysis region [default: 1]
%      N1: end of analysis region [default: length(x)]
%
% Outputs:
% jitter is a structure containing the following values:
%   jitter.PF: the PF of the inter-excitation interval time series
%   jitter.PQ: the PQ of the inter-excitation interval time series
%   jitter.K : the window lengths used to compute the values in jitter.PQ
%
% shimmer is a similar structure, except that the PF and PQ measures are 
% applied to the time series containing the within-excitation (peak-to-peak)
% energy values
%
%        flag:  a Boolean flag set to 1 if no excitations were found 
%               in the analysis region in this case, jitter and shimmer 
%               values are set to 0.
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


if ((nargin < 3) | isempty(Kj)); Kj = [3]; end
if ((nargin < 4) | isempty(Ks)); Ks = [3]; end
x = x(:)';
Kj = Kj(:)';
Ks = Ks(:)';

if ((nargin < 5) | isempty(vexcit)); vexcit = voiced_excite(x,Fs); end
if nargin < 6
  N0=1; N1=length(x);
end
flag = 0;


% Build period_length and energy_within_periods 
% time series to calculate jitter and shimmer respectively
period_length_cum=[]; energy_p2pk_cum=[];
for k=1:vexcit.num_segs
     excite_seg = vexcit.exciteinds{1,k};
     excite_seg = excite_seg(find((excite_seg > N0) & (excite_seg<N1)));	% Retain only indices of interest for the
     period_length_cum = [period_length_cum diff(excite_seg)];			% analysis
     energy=[];
     for cycle=1:length(excite_seg)-1
        energy(cycle) = sum(x(excite_seg(cycle):excite_seg(cycle+1)-1).^2);
     end
     energy_p2pk_cum = [energy_p2pk_cum energy];
end

% Find jitter and shimmer measures (using the perturbation factor and
% perturbation quotient measures for each)
if ~isempty(period_length_cum)
  invalidsize = find (Kj > length(period_length_cum));
  Kj(invalidsize) = 1;  
  [PF_j,PQ_j] = perturb(period_length_cum,Kj);

  invalidsize = find (Ks > length(energy_p2pk_cum));
  Ks(invalidsize) = 1;
  [PF_s,PQ_s] = perturb(energy_p2pk_cum,Ks);
else
  PF_j = 0; PQ_j = zeros(1,length(Kj));
  PF_s = 0; PQ_s = zeros(1,length(Ks));
  flag = 1;
end


jitter.pert_factor = PF_j;
jitter.K = Kj;
jitter.pert_quotnt = PQ_j;
shimmer.pert_factor = PF_s;
shimmer.K = Ks;
shimmer.pert_quotnt = PQ_s; 


