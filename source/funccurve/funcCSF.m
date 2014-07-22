% function [predictedcritvals] = funcCSF(params, xs);
% The contrast sensitivity function.
%
% default: 
% ycrit = p1 + 10.^p3.*(log10(xs) - p2).^2
%
% p1 corresponds to the log contrast threshold at the peak frequency
% p2 corresponds to log10(peak frequency)
% p3 corresponds to log10(inverse bandwidth)
% xs are log10(spatial frequency)
% ycrit is log10(threshold contrast)
% 
% Derivation and description follows:
% The contrast sensitivity function suggested Pelli, Legge and Rubin (JOSA
% 1986) has three parameters: if ycrit is log10(threshold contrast) it implies
%
% ycrit = y0 +K(log10(f) - log10(f0)).^2
% 
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
% 
% Other CSF functions are available (see notes below, and uncomment a line if you prefer that parameterization)
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [predictedcritvals] = funcCSF(params, xs);

% default CSF is a transformation of the Pelle Legge, Rubin, described above

predictedcritvals = params{1} + 10.^params{3}.*(log10(xs) - params{2}).^2;

%%
% Mannos Sakrison, 1974 -- original
% This is the psychometric function from Mannos and Sakrison, 1974 -- 
% !! This returns the sensitivity (1/threshold)!
% CSF(f) = a(b+cf)*exp(-((cf)^d))
% 
% Their estimated parameters were:
% CSF(f) = 2.6*(0.0192 + 0.114*f)*exp(-(0.114*f)^1.1);
% uncomment the line below if you want to use it

% predictedcritvals = params{1}.*(params{2} + params{3}.\xs).*exp(-(params{3}.\xs).^params{4});


%%
% Mannos Sarkison, 1974 -- 3-parameter, log transformed, threshold instead of sensitivity
% This returns log threshold contrast, not sensitivity, for x in cpd.
%
% ycrit = -p1 -log10((10.^p2 + 10.^p3.*xs)) +(10.^p3.*xs)./log(10);
%
% The contrast sensitivity function suggested Mannos & Sakrison, 1974,
% has four parameters and takes the form:
% 
% suitable priors might be centered on 
% log10(260), log10(.0192) log10(.114)
% with all parameters designated for linear
% spacing in the initializing call to fastFull.
% uncomment the line below to use.

% predictedcritvals = -params{1} -log10((10.^params{2} + 10.^params{3}.*xs)) +(10.^params{3}.*xs)./log(10);

%% 
% The contrast sensitivity function of the log-parabola form taken from
% Watson and Ahumada (2005).
% Specifically the function for sensitivity is:
% S(f; f_0, b, a) = {1-a when f<f_0 and S<1-a
%                   {10.^-(log10(f/f_o)./b).^2 otherwise.
% 
% p1 = f_0  the cut-off frequency
% p2 = a    the cut-off sensitivity
% p3 = b    the bandwidth of the CSF.
% 
% uncomment both lines below to use

% predictedcritvals = 10.^-((log10(xs./params{1})./params{3}).^2);
% predictedcritvals(predictedcritvals<(1-params{2}) & xs < params{1}) = (1-params{2});

