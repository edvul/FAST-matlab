% function [predictedcritvals] = funcHyperbolic(params, xs);
%  This is a simple 1-parameter hyperbolic decay (from 1 to 0 for positive
%  x values).  
%   y = 1/(1 + p1*x)
% 
% I can't think of many uses for this -- it was created for a
% delay-discounting experiment.
%
% copyleft Ed Vul and Don MacLeod 2007
%   contact: evul@mit.edu

function [predictedcritvals] = funcHyperbolic(params, xs);

predictedcritvals = 1./(1+params{1}.*xs);