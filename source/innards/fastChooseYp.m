% function [Yp (Yp_sd)] = fastChooseYp(fast, xs, ps);
% Predicts y values for a given p value, (like fastCalcYs
% but fastChooseYp integrates over the whole parameter space.)
% Input:
%   fast structure
%   x an x value
%   P a desired p value.
% Output:
%   Yp: Y value.
%  Optional:
%   Yp_sd: standard deviation of the mean represented in Yp.
% 
% Example:
%   y = fastChooseYp(fast, [1:10], 0.75);
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16


function [varargout] = fastChooseYp(fast, x, p)
    fastSettings;
    if(nargin<3)
        error('Error: fastChooseYp must have both x and p as input.');
    end
    if((numel(x) > 1) && (numel(p)>1) && (numel(x) ~= numel(p)))
        error('Error: x and p must be either single-value or equally-sized vectors.');
    elseif((numel(x) > 1) && (numel(p)>1) && (numel(x) == numel(p)))
        x = reshape(x, 1, length(x));
        p = reshape(p, 1, length(x));
    elseif((numel(x) == 1) && (numel(p) == 1))
        x = x;
        p = p;
    elseif((numel(x) > 1) && (numel(p)==1))
        x = reshape(x, 1, length(x));
        p = repmat(p, 1, length(x));
    elseif((numel(p) > 1) && (numel(x)==1))
        p = reshape(p, 1, length(p));
        x = repmat(x, 1, length(p));
    end
    
    % prior over parameters (P(parameters))
    bigprior = fast.params.core.log10lh - max(fast.params.core.log10lh(:));
    bigprior = 10.^bigprior;
    bigprior = bigprior ./ sum(bigprior(:));
    
    ndp = size(bigprior);
    nd = length(ndp);
    bigprior = repmat(bigprior, [repmat(1,1,nd) numel(x)]);
    
    % p(Y|parameters,P,X)
    ys = fastCalcYs(fast, x, p, fast.params.core.lattice);
    
    % take the expectation of P(Y|P,X) by summing over parameters
    weightedYs = (ys.*bigprior);
    
    muY = sumto1d(weightedYs, nd+1);
    bmuY = repmat(reshape(muY, [repmat(1,1,nd), length(x)]), [ndp, 1]);
    weightedVar = (((ys-bmuY).^2).*bigprior);
    sdY = sqrt(sumto1d(weightedVar, nd+1));
    
    varargout{1} = double(muY);
    if(nargout == 2)
        varargout{2} = double(sdY);
    end