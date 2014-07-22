% function [predictedcritvals] = funcExpS(params, xs);
%  This is a simple 1 parameter exponential decay model... 
% y = exp(-params{1}.*xs);
%
% copyleft Don MacLeod and Ed Vul 2007
%   contact: evul@mit.edu

function [predictedcritvals] = funcExpS(params, xs);

predictedcritvals = exp(-params{1}.*xs);