function downsampled=downsample_to_n_hz(signal,orig_freq,downsampled_freq) 
%
% Usage:
%
% downsampled=downsample_to_n_hz(signal,orig_freq,downsampled_freq)
%
% Low-tech wrapper for resample, making it an error to try to upsample. 
% Downsamples signal to downsampled_freq Hz and returns the resulting 
% signal as a MATLAB vector.  
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


% test sampling rate and downsample to n kHz
if (downsampled_freq>orig_freq) 
     error ('This algorithm assumes an input sampling frequency greater than or equal to n Hz.');
end

ratio = downsampled_freq/orig_freq;
downsampled = resample(signal,round(ratio*1000),1000);
