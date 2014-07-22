% function [varargout] = psyNormal([nchoice],sigmaParameter, yscrit, ys, rs, ps);
% Like all Psychometric functions in fast, it can operate in three modes
% (all psychometric functions used by FAST have to be defined to work in
% these modes):
% 
% [islog] = psyFunction(0);
%   1 (function is defined over (0,Inf)) or 
%   0 (function is defined over (-Inf,Inf)
% [ps] = psyFUNCTION([], parameter, yscrit, ys);
%   output p(response) (no rs or ps provided)
% [ps lls]  = psyFUNCTION([], parameter, yscrit, ys, rs);
%   output p(response) and ll(response) (if rs are provided)
% [ys_p] = psyFUNCTION([], parameter, yscrit, [], [], ps);
%   output y that produces ps(response) (provided) at ycrit (provided)
% 
% Note that here the first parameter is a dummy parameter, in place of 'nchoice' for the NormalN function.
% In calls to this function the first parameter can be empty.
% This fits a complete normal sigmoid, with its inflection point at p = .5,
% to data from YES/NO experiments. (For forced choice, use psyNormal2 or
% psyNormalN.) The stimulus value y could be log or linear--most often it
% could be the signed difference (log or linear) between two stimuli that
% are COMPARED. 
%
% But this function can also be applied to (yes/no) DETECTION experiments.
% If y represents the log of the magnitude of the stimulus to be detected,
% the predicted probability of a positive response goes to zero (leaving
% aside lapses) for zero stimulus magnitude (y - -Inf). (The full sigmoid
% description of detection data is less natural if y is linear with
% stimulus magnitude. There is a rationale for this, if Gaussian noise
% helps a signal linear with y to exceed a high threshold (S.A. Klein,
% 2001), but the predicted false alarm rate at y = 0 can easily be too
% high.) 
%
% In all cases, the returned critical value ycrit is the value associated
% with 50% probability of a positive response,  and sigmaParameter is the
% the standard deviation in y associated with the fitted normal cdf. Thus
% at y = ycrit +sigmaParameter,  ppp = .84.
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [varargout] = psyNormal( nchoice, sigmaParameter, yscrit, ys, rs, ps) % params, xs, ys, rs

AMILOG = 0;

if(AMILOG && (any(yscrit(:)<=0) || any(ys(:)<=0)))
    error('Error: This psychometric function is defined over y=(0,Inf), so it cannot use values where y<=0');
end

if(nargin == 1)
    varargout{1} = AMILOG;
elseif(nargin < 6)
    ppp = 0.5 + 0.5.*erf((ys-yscrit)./(sqrt(2).*sigmaParameter));  % true p, complete normal cdf; sqrt(2) because of erf

    ppp = fastPsyScale(ppp, nchoice, 0);
    varargout{1} = ppp;


    if(nargin == 5)  % we're trying to estimate likelihoods as well as probabilities
        ll = rs.*log(ppp+eps)./2.3026 +(1-rs).*log(1-ppp+eps)./2.3026;
        varargout{2} = ll;
    end
elseif(nargin == 6)
    ps = fastPsyScale(ps, nchoice, 1);

    ys = yscrit + sqrt(2).*sigmaParameter.*erfinv(2*ps-1);
    varargout{1} = ys;
end
