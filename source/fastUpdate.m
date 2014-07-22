% function [fast (resamp)] = fastUpdate(fast, data);
% Updates likelihood lattice, and parameter estimates.
% Input:
%   fast structure
%   data: an Nx3 array of
%       N rows of data, each one contains an
%       [(x value), (y value), (response)]
%   If no data are provided, fastUpdate just updates the parameter
%   estimates.
% Output:
%   fast: updated fast structure.
%   (optional) resamp: 1 or 0: is resampling recommended.
% 
% Example:
%   [fast resamp] = fastUpdate(fast, [10 .6 1]);
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function varargout = fastUpdate(fast, newdata)
    UPDATEPARAMS = 1;
    resamp = 0;
    if(nargin < 2)
        newdata = [];
    end
    if(nargin < 3)
        updatefromscratch = 0;
    end

    if(numel(newdata) > 0) % some new data provided, update with new data
        if(size(newdata, 2) == 3)
            % newdata = newdata;
        elseif(size(newdata,1) ~= 3)
            error('Error: newdata must be Nx3');
        else
            newdata = newdata';
        end
        
        % hard-coded estimate on how many data points are TOO many to
        % update at once -- if too many to update at once, do as a loop
        if((numel(fast.params.core.lattice{1})*size(newdata,1)) > 150000) 
            for i =[1:size(newdata,1)]
                fast = updateData(fast, newdata(i,:));
            end
        else
            fast = updateData(fast, newdata);
        end
        if(UPDATEPARAMS)
            [fast.params.est resamp] = fastCalcEstimates(fast);
        end
    elseif(UPDATEPARAMS) % no new data, must just want to update parameter estimates
        [fast.params.est resamp] = fastCalcEstimates(fast);
    end
    
    % if 'enough' data, suggest resampling, otherwise, don't.
    varargout{1} = fast;
    npp = size(fast.data, 1) / (10*fast.params.n);
    if(nargout == 2)
        if(npp > 3)
            varargout{2} = resamp;
        else
            varargout{2} = 0;
        end
    end
end

%% updateData
function fast=updateData(fast, newdata)
    % update parameter likelihoods
    ycrit = fastCalcYs(fast, newdata(:,1), -1, fast.params.core.lattice);
    
    [ppp dlog10lh] = fastCalcPsLLs(fast, ...
        fast.params.core.lattice{fast.params.n}, ...
        ycrit, ... ycrits
        newdata(:,2),... ys
        newdata(:,3)); % rs
    
    dlog10lh = sum(dlog10lh, fast.params.n+1);
    fast.params.core.log10lh = fast.params.core.log10lh + reshape(dlog10lh, size(fast.params.core.log10lh));
    
    fast.data(end+1:end+size(newdata, 1),:) = newdata;
end
