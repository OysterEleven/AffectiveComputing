function  [loc,val] = findpeak(spec,npicks,rdiff)
%FINDPEAK Finds peaks
% [loc,val] = pickpeak(spec,npicks,rdiff)
%       spec   - data vector or matrix
%       npicks - number of peaks desired              [default = 2]
%       rdiff  - minimum spacing between picked peaks [default = 5]
%       loc    - vector of locations (indices) of the picked peaks
%       val    - vector corresponding values
%       A 0 in location (i,j) of array loc (or a NaN in array val)
%       indicates that the j-th data vector has less than i peaks
%       with a separation of rdiff or more.

% ---- parameter checks  -------------------------------------------

if (exist('rdiff') ~= 1)  rdiff =  5;                  end
if (exist('npicks') ~= 1) npicks = 2;                  end

% ---- convert row vectors to col vectors  -------------------------

[mrows,ncols]  = size(spec);
if (mrows==1) mrows=ncols; ncols=1; spec = spec(:);   end

% ---- edit out NaNs and Infs ---------------------------------------

good = find (finite(spec));
rmin = min(spec(good)) - 1;
bad  = find(~finite(spec));
if (~isempty(bad))
   spec(bad) = ones(size(bad)) * rmin;
end

% ---- find a peak, zero out the data around the peak, and repeat

val =  ones(npicks,ncols) * NaN ;
loc =  zeros(npicks,ncols) ;

for k=1:ncols
                                           % Find all local peaks:
    dx = diff([rmin; spec(:,k); rmin]);    % for a local peak at either end
    % The following line of code is buggy, as it ANDS two conditions that
    % include the case of the slope being strictly zero at both edges of a
    % point being considered ==> This is not a peak; this is a flat line!
    % 
    % Instead, consider 3 cases for a peak separately
    % i. strictly +ve slope followed by strictly -ve slope
    % ii. strictly +ve slope followed by zero slope
    % iii. zero slope followed by strictly -ve slope 

    % (This line of code was included in the Matlab routine, and it's commented out
    % lp = find(dx(1:mrows)   >= 0 ...
    %        & dx(2:mrows+1) <=0);          % peak locations

    dx_sgn = sign(dx);
    lp_plus_minus = find(dx_sgn(1:mrows)   == 1 ...
                      & dx_sgn(2:mrows+1) == -1);
    lp_plus_zero = find(dx_sgn(1:mrows)   == 1 ...
                      & dx_sgn(2:mrows+1) == 0);
    lp_zero_minus = find(dx_sgn(1:mrows)   == 0 ...
                      & dx_sgn(2:mrows+1) == -1);
    lp_all = unique([lp_plus_minus; lp_plus_zero; lp_zero_minus]);

    vp = spec(lp_all,k);                       % peak values

    for p=1:npicks
       [v,l] = max(vp);                   % find current maximum
       val(p,k) = v;  loc(p,k) = lp_all(l);   % save value and location

       ind = find(abs(lp_all(l)-lp_all) > rdiff);  % find peaks which are far enough away

       if (isempty(ind))
           break                           % no more local peaks to pick
       end
       vp  = vp(ind);                      % shrink peak value array
       lp_all  = lp_all(ind);                      % shrink peak location array
    end
end

