% function leftsum=sumto1d(leftsum, d)
% you would think this functionality would be included in Matlab, but no.
% It takes an N dimensional array (leftsum), and sums it to 1 dimension, as
% specified by d.
% 
% Used for a number of purposes -- e.g., marginalization
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function leftsum=sumto1d(leftsum, d)
    for i=[1:length(size(leftsum))]
        if(i~=d)
            leftsum = sum(leftsum,i);
        end
    end
    leftsum = squeeze(leftsum);