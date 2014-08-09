% function Y = fastChooseY(fast, xs, method);
% A wrapper function for the basic function of:
%   fastChooseYp, fastChooseYent, fastChooseYpre, fastCalcYs
% Use these functions directly for speed, and to make them more efficient:
% If you specify their parameters yourself, you may achieve greater
% efficiency from testing than the defaults we try to figure out here.
% Also, getting to know and using those functions directly allows you to
% use a few additional, useful features on top of simply choosing the next
% stimulus.
% 
% Input:
%   fast structure
%   xs: 1 or a vector of x values
%   method: which way of choosing a Y value? 
%       (default)   'int': integrate over all parameter values.
%                       fastChooseYp
%                   'ent': choose expected minimum entropy setting.
%                       fastChooseYent
%                   'pre': expected minimum predictive variance setting.
%                       fastChooseYpre
%                   'mle': use best point parameter estimate.
%                       fastCalcYs
% Output:
%   Y: Y value.
% 
% Example:
%   y = fastChooseY(fast, 4);
% 
% A note on defaults, fastChooseY tries to come up with a sensible Y range
% or P value for you to be sampling, but you may be wiser to use the
% functions directly and specify these parameters yourself.
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function Y = fastChooseY(fast, x, method)
    if(nargin < 3)
        method = 'int';
    end
    
    switch(method)
        case {'int'} % fastChooseYp(fast, x, p)
            % first, which p?
 			ps = [.4 .5 .6];
            ps = fastPsyScale(ps(ceil(rand(1,1)*length(ps))), fast.params.nchoice, 0); % scale as needed
            
            Y = fastChooseYp(fast, x, ps);
            
        case {'mle'} % fastCalcYs(fast, x, p)
            % first, which p?
			ps = [.4 .5 .6];
            ps = fastPsyScale(ps(ceil(rand(1,1)*length(ps))), fast.params.nchoice, 0); % scale as needed
            
            Y = fastCalcYs(fast, x, ps);
            
        case {'ent'} % fastChooseYent(fast, x, y, iterations)
            % first, which Y range?
            ps = [.1 .9]; % take ys for 10 to 90% scaled
            ps = fastPsyScale(ps, fast.params.nchoice, 0); % scale as needed
            ys = fastCalcYs(fast, repmat(x, length(ps), 1), ps);
            
            Y = fastChooseYent(fast, x, ys);
            
        case {'var'} % fastChooseYpre(fast, x, y, iterations)
            % first, which Y range?
            ps = [.1 .9]; % take ys for 10 to 90% scaled
            ps = fastPsyScale(ps, fast.params.nchoice, 0); % scale as needed
            ys = fastCalcYs(fast, repmat(x, length(ps), 1), ps);
            
            Y = fastChooseYpre(fast, x, ys);
    end
    