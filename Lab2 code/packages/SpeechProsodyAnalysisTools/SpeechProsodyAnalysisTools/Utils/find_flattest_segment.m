function [first, last, seed]=find_flattest_segment(signal, k, left_mult, right_mult)
% Usage: 
%
% [first, last, seed]=find_flattest_segment(signal, k, left_mult, right_mult);
%  
% First finds the "flattest" length-k segment of the signal, that is, the  
% segment for which the sum of the absolute value of the first difference
% is minimum (seed is the index of the beginning of this length-k segment), 
% then grows this region to the right and left using statistical
% criteria based on the mean and standard deviation of the growing segment. 
% 
% In particular, the region is grown sample by sample to the right as long
% as each new sample lies within right_mult standard deviations of the mean of 
% the growing region. (That is, the mean and standard deviation are recalculated 
% every time a new sample is admitted to the region.)
%
% The region is then grown to the left sample by sample using a slightly different 
% criterion. Each new sample must still lie within left_mult standard deviations of 
% a mean, but the mean and standard deviation used are the final mean and 
% standard deviation resulting from the rightward expansion; they are not 
% recalculated following the admittance of new samples.
%
% First and last are, appropriately enough, the first and last index of the 
% final region.
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


if nargin==1 
  k=5;
  left_mult=2;
  right_mult=2;
end
if nargin==2
left_mult=2;
right_mult=2;
end
if nargin==3
right_mult=2;
end
len=length(signal);
%
% Calculate change function.
%
if len<k; k = len; end
for i=1:(length(signal)-(k-1))
change_function(i)=sum(abs(diff(signal(i:i+(k-1))))); 
end
[minval,index]=min(change_function);
%
% Initial length-k region.
%
seed=index(1); % if >1 minimal elements, just grab the first
first=seed;
last=first+(k-1);
%
% Grow region to the right.
%
mn=0; % inconsequential value for variable persistence
stddev=0; % ditto
moreToDo=1;
while((last+1<len) & ... % avoid running off the end 
      moreToDo==1)       % |(mean-point)|<right_mult*stddev
mn=mean(signal(first:last));
stddev=std(signal(first:last));
next_sample_right=signal(last+1);
if ( abs(next_sample_right - mn) < (right_mult*stddev) )
     last=last+1;
else moreToDo=0;
end%if/else
end%while
%
%fprintf(['Standard deviation is ' num2str(stddev) '.\n.']);
%
% Now grow to the left. Use final mean and stddev from rightward 
% expansion, which remain in mn and stddev.
%
while((first>1) & ...  % avoid running off the front
      abs(signal(first)-mn)<(left_mult*stddev))
first=first-1;
end%while
%
% Make a pretty picture
%
%figure;
%plot(signal);
%seedvec=zeros(1,len);
%seedvec(index)=.25;
%outputvec=zeros(1,len);
%outputvec(first:last)=.15;
%hold on
%stem(outputvec,'r')
%     stem(seedvec,'g')
%zoom on
