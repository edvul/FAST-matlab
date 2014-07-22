% function [varargout] = funcCalcPsLLs(fast, parameter, yscrit, ys, rs); 
%   This is an all purpose function for interacting with the psychometric
%   function predictions.
%   For the most part, this should not be used directly -- it is primarily
%   for internal FAST functions.  However, usage description follows.
%
%   Currently it operates under three schemes, which are determined by the
%   type and arrangement of parameters provided.
% 
%   Operating scheme 1: Produce probabilities of responses
%    Input: 
%     fast: the parametric staircase structure
%     PARAMETER: the psychometric function parameter for this calculation
%     YSCRIT: a vector of critical y values
%     YS: a vector of y values (corresponding to the critical values above)
%    Output:
%     ppp: probability of a response for each of the ycrit,y pairs.
% 
%   Operating scheme 2: Produce probabilities and likelihoods for arrays 
%       of ycrit, ys, rs
%    Input:
%     fast: the parametric staircase structure
%     PARAMETER: a single psychometric function parameter
%     YSCRIT: an array of critical y values
%     YS: an array of y values(corresponding to the critical values above)
%     RS: an array of responses for the above ys, yscrit.
%    Output:
%     ppp: probability of a response for each of the ycrit,y sets.
%     lls: log likelihood of each ycrit,y,r set
%    Optional, but not recommended (only used within the FAST code, and
%    is not be guaranteed to do what you predict): you can put in an array
%    of yscrit and only a vector of ys, and rs;  these are then filled to
%    size with yscrit (but, there are assumptions made here about how these
%    are all sized, and so your best bet would be to properly size all of
%    the parameters yourself).
% 
%   Operating scheme 3: Produce probabilities and likelihoods for a large set of 
%       psychometric function parameters.
%    Input:
%     fast: the parametric staircase structure
%     PARAMETER: a vector of psychometric function parameters for
%       comparison
%     YSCRIT: a vector of critical y values
%     YS: a vector of y values (corresponding to the critical values above)
%     RS: a vector of responses for the above ys, yscrit.
%    Output:
%     ppp: probability of a response for each of the ycrit,y,psychparam sets.
%     lls: log likelihoods of each ycrit,y,psychparam set
%     
%
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16


function [varargout] = fastCalcPsLLs(fast, parameter, yscrit, ys, rs) 

if(nargin == 4) % we are just polling for probabilities for ycrit,y pairs, parameter should be fixed.
    if(length(parameter) > 1)
        error('Error: funcLogistic() either takes 4 arguments and a fixed parameter, or 5 arguments and variable parameters');
    end
    
    % this only works with one parameter;
    parameter = squeeze(parameter);
%     bigparam = repmat(parameter, size(yscrit));
    bigyscrit = yscrit;
    
    bigys = reshape(ys, size(yscrit));

    ppp = fast.func.psych(fast.params.nchoice, parameter, bigyscrit, bigys);
    varargout{1} = ppp;
end

if(nargin == 5)  % we're trying to estimate likelihoods as well as probabilities
    ndp = size(yscrit); % figure out if this should return likelihood of psych parameters, or of curve params with fixed psych params
    np = fast.params.n;
    if(np < length(ndp))
        bigparam = repmat(parameter, [repmat(1, 1, np) ndp(end)]);
    else
        bigparam = parameter;
    end
    bigyscrit = yscrit;
    if(numel(ys) == 1) % only 1 y... this here to remedy problems that arise with 1d matrices.
        nd = length(ndp);
%         disp a
        bigys = repmat(reshape(ys, [repmat(1, 1, nd), length(ys)]), [ndp(1:nd), 1]);
        bigrs = repmat(reshape(rs, [repmat(1, 1, nd), length(rs)]), [ndp(1:nd), 1]);
    elseif(numel(ys)==numel(yscrit)) %do no resizing, everything should be perfectly presized;
%         disp b
        bigys = ys;
        bigrs = rs;
    elseif(numel(ys) == length(ys)) % ys are a vector assume rs are too and that yscrit are coming from parameter estimation
%         disp c
        nd = length(ndp) - 1;
        bigys = repmat(reshape(ys, [repmat(1, 1, nd), length(ys)]), [ndp(1:nd), 1]);
        bigrs = repmat(reshape(rs, [repmat(1, 1, nd), length(rs)]), [ndp(1:nd), 1]);
    end
    [ppp ll] = fast.func.psych(fast.params.nchoice, bigparam, bigyscrit, bigys, bigrs);
    varargout{1} = ppp;
    varargout{2} = ll;
end