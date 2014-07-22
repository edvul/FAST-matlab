% function [predictedcritvals] = funcExp(params, xs);
%  an exponential decay/rise function
% y = p1 + (p2-p1)*(1-exp(-x/p3));
% 
% Parameters may be easily interpreted as:
% p1: initial value (when x=0)
% p2: asymptote approached as x >> Inf
% p3: scaling (usually, time) constant
% 
% copyleft Ed Vul and Don MacLeod 2007
%   contact: evul@mit.edu

function [predictedcritvals] = funcExp(params, xs);

predictedcritvals = params{1} + (params{2} - params{1}).*(1-exp(-xs./params{3})); 