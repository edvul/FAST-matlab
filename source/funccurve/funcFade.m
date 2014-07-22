% function [predictedcritvals] = funcFade(params, xs);
%  to vary 2 parameters of a decay multiplier, ranging 1 to 0
% In addition to the initial and final values (1 and 0), the 1/e point is pinned; 
% Curve is exp(-f(xs)) where f(xs) is different for Versions 1 and 2. 
% Version 1: Curve is exp(-f(x)) where f(x)= px + (1-p).x^2
% p = 0 makes it Gaussian, but for p > 0, < 1 there is a linear initial descent
% Rate of descent is proportional to p for a given time constant
% For p > 1 f(x) is not monotonic, and the 'decay' is reversed at long times.
% p = 1 gives an exponential, bigher values are needed to mimic sum of exps.
% p can be as high as 1.2 and give monotonic decay for about 3 time constants.
%
% Version 2 has a better focus on heavy-tailed decay: f(x) = log(1+x)/log(2), which doesn't
% reverse, for 1 > p > 0; exp for p = 1, heavy-tailed distortion of exponential for p < 1,
% Range could be from -.3 or so (monotonic decay to 5 time constants or so) 
% up to high values (heavy peaked but still monotonic at high values). 
% The upper limit for monotonicity at t = 0 is 
% (1-p)/log(2) = -p, or -p = 1.4427(1-p), or .4427p = 1.4427, or p =  3.259
% Lower limit is strictly p = 0, but values down to -.3 or so could be useful..


% copyleft Don MacLeod and Ed Vul 2007
%   contact: evul@mit.edu
function [predictedcritvals] = funcFade(params, xs);
if iscell(params)
    params = cell2mat(params);
end
p = params(1);
t = xs/params(2);
ft = p.*t +(1-p).*log(1+t)./log(2);  % version (2)
% Version 1 ft = p.*t +(1-p).*t.^2
predictedcritvals = exp(-ft);
