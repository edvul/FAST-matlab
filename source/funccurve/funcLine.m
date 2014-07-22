% function [predictedcritvals] = funcLine(params, xs);
% A line.
% predictedcritvals = params{1} + params{2}.*xs;
%
% copyleft Ed Vul & Don MacLeod, 2007
% contact: evul@mit.edu
% version: FAST v2.5

function [predictedcritvals] = funcLine(params, xs);

predictedcritvals = params{1} + params{2}.*xs;