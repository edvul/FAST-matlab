% function [predictedcritvals] = funcCSF(params, xs);
% The contrast sensitivity function of the log-parabola form taken from
% Watson and Ahumada (2005).
% Specifically the function for sensitivity is:
% S(f; f_0, b, a) = {1-a when f<f_0 and S<1-a
%                   {10.^-(log10(f/f_o)./b).^2 otherwise.
% 
% p1 = f_0  the cut-off frequency
% p2 = a    the cut-off sensitivity
% p3 = b    the bandwidth if the CSF.


function [predictedcritvals] = funcCSFs(params, xs);

predictedcritvals = 10.^-((log10(xs./params{1})./params{3}).^2);
idx = predictedcritvals<(1-params{2}) & xs < params{1};
predictedcritvals(idx) = (1-params{2}(idx));
                
                