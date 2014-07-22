% function [predictedcritvals] = funcExp(params, xs);
% An unlimited-degree generalized polynomial function.  
% 
% The order of the polynomial is determined by the number of parameters.
%   1 parameter is a 0th order polynomial (a singular value):
%       y = p1
%   2 parameters is a 1st order polynomial (a line):
%       y = p1 + p2*x
%   3 degrees is a 2nd order polynomial (quadratic):
%       y = p1 + p2*x + p3*x^2
%   and so forth; n parameters produces an (n-1)th order polynomial:
%       y = SUM_n{ pn * x^(n-1) }
% 
% copyleft Ed Vul and Don MacLeod 2007
%   contact: evul@mit.edu


function [predictedcritvals] = funcPolynomial(params, xs);

predictedcritvals = zeros(size(xs));
for i=[1:length(params)]
    predictedcritvals = predictedcritvals + params{i}.*(xs.^(i-1));
end