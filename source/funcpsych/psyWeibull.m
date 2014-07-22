% function [varargout] = psyWeibull(nchoice, parameter, yscrit, ys, rs, ps);
% 
% Basic form: (1-exp(-(Y/Y*)^B))
% Use: Weibull functions are defined over the range [0:+Inf).
%   Thus, the Y dimension is assumed ot be logarithmically scaled 
%   (e.g., contrast) which means that Y values can NOT be negative.
%   Since the function operates over Y/Ycrit (rather than Y-Ycrit) 
%   like most other psychometric functions, it is translation invariant 
%   in log-space, which is a desirable feature for many psychophysical 
%   measures.
%   The slope parameter is commonly referred to as 'Beta' and usually is
%   somewhere between 2 and 3.5.
% 
% Like all Psychometric functions in FAST, it can operate in three modes
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
% All psychometric functions in FAST are defined as full sigmoidal
% functions spanning the range of [0 1].  They are then scaled based on the
% experimental paradigm (nAFC, matching, detection), to fit the relevant
% range of probabilities (e.g., .5 to 1 for 2 AFC) using the function
% fastPsyCrop based on the nchoice parameter and some preset lapse
% parameters.
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [varargout] = psyWeibull(nchoice, parameter, yscrit, ys, rs, ps); % params, xs, ys, rs

AMILOG = 1;

if(nargin > 1 && AMILOG && (any(yscrit(:)<=0) || any(ys(:)<=0)))
    error('Error: This psychometric function is defined over y=(0,Inf), so it cannot use values where y<=0');
end


if(nargin == 1)
    varargout{1} = AMILOG;
elseif(nargin < 6)
    ppp = 1-exp(-(ys./yscrit).^parameter);  % usee /max(yscrit,eps) to avoid NaN??
    ppp = fastPsyScale(ppp, nchoice, 0);
    varargout{1} = ppp;

    if(nargin == 5)  % we're trying to estimate likelihoods as well as probabilities
        ll = rs.*log(ppp+eps)./2.3026 +(1-rs).*log(1-ppp+eps)./2.3026;
        varargout{2} = ll;
    end
elseif(nargin == 6)
    ps = fastPsyScale(ps, nchoice, 1);
    
    ys = yscrit.*(-log(1-ps)).^(1./max(parameter,eps));
    varargout{1} = ys;
end
 