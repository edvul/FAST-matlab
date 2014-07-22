% function [predictedcritvals] = funcMscale(params, xs);
% This is really a line, but it is re-parametrized so that the parameters
% have a more intuitive interpretation as Magnification scaling (with
% eccentricity)
%   y = p1*(1 + p2*x)
% 
% thus, p1 is the foveal magnificationm, and p2 corresponds to the
% proportional increase per unit of eccentricity (x)
%
% copyleft Don MacLeod and Ed Vul 2007
%   contact: evul@mit.edu

function [predictedcritvals] = funcMscale(params, xs);

predictedcritvals = params{1}.*(1 + params{2}.*xs);