% function [predictedcritvals] = funcCSF(params, xs);
% This returns log threshold contrast, not sensitivity, for x in cpd.
% The contrast sensitivity function suggested Mannos & Sakrison, 1974,
% has four parameters and takes the form:
% 
% CSF(f) = a(b+cf)*exp(-((cf)^d))
% 
% For contrast sensitivity expressed as the reciprocal of threshold percent
% contrast, and f in cycles per degree, their estimated parameters were:
% CSF(f) = 2.6*(0.0192 + 0.114*f).*exp(-(0.114*f).^1.1);
% For contrast from 0 to 1, the value of a corresponding to
% Mannos and Sakrison's 2.6 should be 2.6*100 or 260.
% This gives a lowest contrast threshold of about .01 at about 8 cpd.
% Generally the psychometric function will be approximately translation-
% invariant on a log scale, so critical values for FAST should be given as
% log10(contrast). If the Mannos and Sakrison formula is adopted, the 
% critical value is then as below, and for Mannos and Sakrison's parameter
% values, the FAST parameter vector is params = [ 260  .0192  .114  1.1].
% But the first two parameters could reasonably be set to be handled as 
% log values in FAST so that lattices have equal spacing in the logs;
% their priors would then be specified in terms of their
% logs in initializing the FAST structure. This is done in funcCSF, which
% accepts parameters that are exponentiated from the values used within
% FAST for parameter lattice creation and averaging.
% 
% Easier: make all parameters linear ones in the sense that they are 
% not only treated linearly in FAST but are passed directly to funcCSF. 
% This entails writing the function in the final form below. In this case
% suitable priors might be centered on log10(260), log10(.0192),
% log10(.114) and log10(1.1), with all parameters designated for linear
% spacing in the initializing call to fastFull.
% This script, funCSFlogMS, does that. Probably this should be the 
% default CSF function, (or the PelliLeggeRubin function could be good).
% Actual parameter values vary greatly with conditions.
% parameter 1 is proportional to peak contrast sensitivity (but roughly e 
% times greater).
% parameter 2 sets the zero-frequency sensitivity as a
% fraction of the peak.
% parameter 3 is very roughly the reciprocal of the best frequency.
% parameter 4 is 1 for exponential asymptotic high-frequency CSF (a 
% straight line in a plot of log(cs) vs f), and 2 for a Gaussian asymptote
% (a parabola in a semilog plot). Should be close to 1.

% copyleft Ed Vul and Don MacLeod 2007
%   contact: evul@mit.edu

function [predictedcritvals] = funcCSFlogMS3(params, xs);
% old format, with FAST's log values exponentiated in passing to func:
% predictedcritvals = -log10(params{1}) -log10((params{2} + params{3}.*xs))...
%    +(params{3}.*xs)./log(10); % funcCSF.m

% or better, if the 1st three  parameter values passed by fast are the
% decimal logs of the quantities used above (so that in this case all 
% passed parameters can be treated as linear ones within FAST and are 
% passed directly to funcCSFlogMS): 
predictedcritvals = -params{1} -log10((10.^params{2} + 10.^params{3}.*xs))...
   +(10.^params{3}.*xs)./log(10);