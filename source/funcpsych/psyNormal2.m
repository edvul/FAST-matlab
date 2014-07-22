% function [varargout] = psyNormal2(nchoice, parameter, yscrit, ys, rs, ps);
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
% Uses: Normal2 applies where Gaussian noise is assumed in (only) a 2AFC
% task. When Gaussian noise contaminates a signal linear with stimulus
% strength, the implied psychometric function is the upper half (only) of a
% normal cdf, hence not a sigmoid. The stimulus strength is scaled by the
% critical value, here defined as the value yielding 84% correct, a
% standard normal deviate of +1 for the difference in signals between the
% two intervals. Empirical 2AFC psychometric functions, however, are
% sigmoidal: d', the effective normal deviate, increases as a power
% function of stimulus strength with exponent between 2 and 3 or so
% (Nachmias and Kocher, 1974). This exponent can serve as a second variable
% parameter (in addition to ycrit), that determines the slope of the
% psychometric function. When probability of a positive response (ppp) is
% expressed as a function of log(y), the slope at a given ppp is
% independent of ycrit, since the curves differ only by a translation in
% log(y), but is proportional to the exponent parameter. The ppp is defined
% over the range of linear stimulus strengths [0:+Inf), thus it won't work
% if Y values can be negative. The slope parameter is commonly referred to
% as 'Beta' and usually is somewhere between 2 and 3.5. Note that although
% p increases as a power of stimulus strength for this function and the
% Weibull function, and also for the Logistic function fitted to a log
% stimulus strength abscissa, the shapes of these functions differ
% considerably as they approach the high p asymptote, the logistic function
% being the shallowest (for a given initial exponent) and the normal2
% function the steepest. 
%
% For consistency with usage for other psychometric functions,
% nchoice is the first parameter, but it must either be = 2 or empty.
%
% Alternatively, the COMPLETE normal cdf could be applied to NAFC data by
% adopting the LOG of the stimulus strength as the independent variable
% (y), with a guessing correction to locate the ppp at 1/nchoice at y =
% -Inf. This can give good fits if position and slope of the normal cdf are
% both free parameters. But the rationale in terms of Gaussian noise is
% lacking. This is what is done in psyNormalN. 
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [varargout] = psyNormal2( nchoice, parameter, yscrit, ys, rs, ps); % params, xs, ys, rs


AMILOG = 1;

if(AMILOG && (any(yscrit(:)<=0) || any(ys(:)<=0)))
    error('Error: This psychometric function is defined over y=(0,Inf), so it cannot use values where y<=0');
end
    
lapseprob = .03; % .985 asymptote
if(nchoice ~=2)
    disp('Error: 1st parameter, nchoice is a placeholder in psyNormal2, can only be 2 or []');
end

nchoice = 2; % nchoice is a dummy parameter ## fix

if(nargin == 1)
    varargout{1} = AMILOG;
elseif(nargin < 6)
    ppp = 0.5 + 0.5.*erf(((ys./yscrit).^parameter)./sqrt(2));  % compressed for guessing, p = .5 to 1   for y >=0; sqrt(2) because of erf
    ppp = (1-lapseprob).*ppp +lapseprob/2;  % lapse correction: was as below, rewritten for clarity
%    ppp = 1./nchoice + (1-lapseprob).*(ppp - 1./nchoice);

      varargout{1} = ppp;


    if(nargin == 5)  % we're trying to estimate likelihoods as well as probabilities
        ll = rs.*log(ppp+eps)./2.3026 +(1-rs).*log(1-ppp+eps)./2.3026;
        varargout{2} = ll;
    end
elseif(nargin == 6) % probability of pos. response given, want corresponding ys:
    ps = (ps-(1./nchoice))./(1-lapseprob)+(1./nchoice); % should we have this here?
    if ps == 0
        ys = -Inf;
    elseif ps ==1
        ys = Inf;
    else
% ((y/yscrit).^parameter)./sqrt(2) =  erfinv(2*(p - .5))
% (y/yscrit).^parameter = sqrt(2)*erfinv(2*(p - .5))
% (y/yscrit) =   (sqrt(2)*erfinv(2*(p - .5))).^(1/parameter)  % ******* NB sqrt is also  to 1/param
% y =   yscrit.*(sqrt(2)*erfinv(2*(p - .5))).^(1/parameter) 

    ys = (sqrt(2)*erfinv(2*ps - 1)).^(1./parameter).*yscrit;  
    end
    varargout{1} = ys;
end
