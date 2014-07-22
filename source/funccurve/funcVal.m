% function [predictedcritvals] = funcVal(params, xs);
%  This is a simple 1 parameter model... makes FAST into a one-dimensional 
%   staircase
%  (although one that updates the slope parameter and has other FAST 
%  functionality)
%
% copyleft Ed Vul & Don MacLeod, 2007
% contact: evul@mit.edu
% version: FAST v2.5

function [predictedcritvals] = funcVal(params, xs);

predictedcritvals = params{1}.*(xs.^0);