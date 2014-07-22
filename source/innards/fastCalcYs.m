% function [ycrit] = fastCalcYs(fast, xs, ps, params);
% Predicts y values for a given x and p value.
% Input:
%   fast structure
%   XS a vector or xs
%   Ps: either one p, or Ps the same size as the xs vector:
%       this indicates which y values to generate (i.e. 
%       those that would produce a particular p(response)
%       rather than the usual default critical value.
%       if set to -1 it just produces the 'critical value' used by the
%       psychometric function...
% Optional:
%   PARAMS can be one of two things:
%       a cell array of parameters
%    or a string indicating which stored parameter set to use:
%       'latticeMax': best with little data
%       'margMean': totally useless, as it turns out
%       'margXCMean': totally useless, as it turns out
%       'gaussMean': Not yet fully functional ##
%       (WITHOUT THIS IT DEFAULTS TO gaussMean)
% 
% Example:
%   y = fastPredict(fast, [1:10], 0.75);
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [ycrit] = fastCalcYs(fast, xs, ps, params)
if(nargin < 3)
    error('fastCalcYs needs three parameters: fast structure, X value, desired P value');
end

if(nargin == 4)
    if(iscell(params))
        params = params;
    elseif(~isempty(params))
        switch(params)
            case {'latticeMax'}
                params = fast.params.est.('latticeMax');
            case {'gaussMean'}
                params = fast.params.est.('gauss').('mean');
            case {'margMean'}
                params = fast.params.est.('marg').('mean');
            case {'margXCMean'}
                params = fast.params.est.('margXC').('mean');
            otherwise
                error(sprintf('Provided parameter set (%s) not recognized', params));
        end
    else
        params = fast.params.est.('gauss').('mean');
    end
else
    params = fast.params.est.('gauss').('mean');
end
    
nel = numel(params{1});
ndp = size(params{1});
nd = fast.params.n;
for i=[1:length(params)]
    if(nel == 1)
        bigparams{i} = repmat(params{i}, 1, length(xs));
    else
        bigparams{i} = repmat(params{i}, [repmat(1, 1, nd), length(xs)]);
    end
end
clear tmp;
if(nel > 1)
    tmp = reshape(xs, [repmat(1, 1, nd), length(xs)]);
    bigxs = repmat(tmp, [ndp, 1]);
else
    bigxs = reshape(xs, 1, length(xs));
end

ycrit = fast.func.curve(bigparams(1:fast.params.n-1), bigxs);
if(ps ~= -1)
    if(numel(ps) == 1)
        ps = repmat(ps, size(ycrit));
    elseif(numel(ps) ~= numel(xs))
        error('ps must be same size as xs, or just one value');
    else % reshape to fit number of parameters
        if(nel >1)
            ps = reshape(ps, [repmat(1, 1, nd), length(xs)]);
            ps = repmat(ps, [ndp, 1]);
        else
            ps = reshape(ps, 1, length(xs));
        end
    end
    ycrit = fast.func.psych(fast.params.nchoice, ...
                         bigparams{fast.params.n}, ...
                         ycrit, ...
                         [], ...
                         [], ...
                         ps);
    ycrit = double(ycrit);
end