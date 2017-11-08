function B=stft(signal, windowvector)
%%%%%%%%%%%%%%%%%%%%% STFT.M %%%%%%%%%%%%%%%%%%%%%%%%%
%
% Special-purpose implementation of the short-time 
% Fourier transform, designed primarily for use with 
% excitation_instants.m. An STFT analysis window of 
% size w=length(windowvector) is centered on every 
% sample of the input signal. That is, successive
% STFT frames overlap by the maximum sensible number 
% of samples, and one transform is calculated for each 
% sample of the input signal. The signal is padded in 
% front and in back with (w-1)/2 zeros. Also, the length 
% of the FFT calculated for windowed segments is always 
% the lowest power of two exceeding twice the segment 
% length. 
%
% For a more general STFT implementation, see specgram.
% 
% USAGE
%
% B=stft(signal, windowvector);
%
% where:
% -B is a (complex) matrix of size fftlength x signallength
%  whose columns contain the Fourier transforms of windowed
%  signal segments;
% 
% -signal is the signal to be STFTed, and 
%
% -windowvector is a vector of windowing coefficents.
%  Use the window('type', length) command to create such
%  a vector. The length of windowvector implicitly defines
%  the frame size and fftlength used in the calculation and
%  *MUST BE ODD*.
%
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



signal=signal(:)'; %' signal is a ROW vector  
windowvector=windowvector(:); % make the windowvector a COLUMN
%
% Length calculations
%
signallength=length(signal); % length of UNPADDED SIGNAL
windowlength=length(windowvector);
if (mod(windowlength,2)==0)
     error('STFT: windowvector must have ODD length!');
end%if
fftlength=2^ceil(log2(windowlength*2));
%
% Screen output
%
fprintf(['STFT> signal length: ' num2str(signallength) ' samples.\n']);
fprintf(['STFT> window size: ' num2str(windowlength) ' samples.\n']);
fprintf(['STFT> length of FFT: ' num2str(fftlength) ' samples.\n']);
%
% Zero-pad input signal
%
paddedsignal=[zeros(1,(windowlength-1)/2) signal zeros(1,(windowlength-1)/2)];
%
% Indices for voodoo magic indexing 
%
rowindex=[1:windowlength]'; %' row index is a COLUMN
colindex=1:signallength; % use UNPADDED signal length
%
% Aforementioned voodoo magic to create matrix of signal segments
%
signalmatrix=zeros(windowlength,signallength); 
signalmatrix(:)=paddedsignal(rowindex(:,ones(1,signallength))+colindex(ones(windowlength,1),:)-1); % Damn Gina! Props to the specgram guy!
%
% Apply the window to each segment of the padded signal
%
windowedsignalmatrix=windowvector(:,ones(1,signallength)).*signalmatrix;
%
% Calculate the FFT of each windowed segment
%
% (here, FFT will automatically zero-pad the whole matrix such 
% that each windowed segment of the signal is at least twice 
% the frame [window] size) 
%
B=fft(windowedsignalmatrix,fftlength,1);
%





