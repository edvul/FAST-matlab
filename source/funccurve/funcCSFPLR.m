% function [predictedcritvals] = funcCSFPLR(params, xs);
% The contrast sensitivity function suggested Pelli, Legge and Rubin (JOSA
% 1986) has three parameters: if ycrit is log10(threshold contrast) it implies
% ycrit = y0 +K(log10(f) - log10(f0)).^2
% where y0 is the log contrast threshold at the peak frequency f0, and K is
% an inverse bandwidth parameter. This is not meant to fit the CSF all the
% way to the resolution limit.
% 
% Generally the psychometric function will be approximately translation-
% invariant on a log contrast scale, so critical values returned to FAST 
% by funcxy functions like this one should be given as log10(contrast). 
%
% It would be natural to sample at linear intervals in y0, log10(f0) and
% log10(K). To do this using y0, f0 and K as the passed values,
% f0 and K could be set to be handled as  log values in FAST, 
% with priors specified in terms of their logs when initializing
% the FAST structure.
% Here we take an alternative approach, passing y0, log10(f0) and log10(K)
% as arguments (recall that y0 is log threshold contrast).
% 
% Actual parameter values vary greatly with conditions.
% Parameter 1 is y0, the peak of ycrit = log10(contrast threshold). 
% Can be treated in FAST as linear, with prior specified directly in
% log10(contrast). Thus if using raw (not %) contrast, the prior for y0
% might be centered near -2.
% Parameter 2 is the corresponding log10(f0),also treated as linear.
% Parameter 3 is log10(K), with K = 1/bandwidth in log10(frequency); 
% log10(K)can be negative, so can be treated by FAST as a linear parameter. 
% So ycrit can be written as y0 +10^log10(K)*(log10(f)-log10(f0))^2
% or params{1} + 10.^params{3}.*(log10(f) - params{2}).^2
% xs could be f or log10(f), but the choice below is xs = f (vector in cpd)
% 
% If y0 = -2, f0 = 8, K = 5 for 1% threshold at 8 cpd, 
% params = {-2  0.9 0.7}
% corresponding cell array of agnostic priors for FastFull
% (add psychometric parameter): 
% {[0 -2 1] [0 .9 1] [ 0 0.7 1]}

% copyleft Ed Vul and Don MacLeod 2007
%   contact: evul@mit.edu

function [predictedcritvals] = funcCSF(params, xs);

predictedcritvals = params{1} + 10.^params{3}.*(log10(xs) - params{2}).^2;
