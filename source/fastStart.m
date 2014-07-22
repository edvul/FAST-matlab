% (For an in-depth description of FAST toolbox, please see the manual)
%
% function fast = fastStart(nchoice, funccurve, parameters)
%
% fastStart() is a simplified wrapper function for fastFull().
% Both initialize a FAST structure; however fastStart is simpler.
% Input:
%   nchoice: [int] task (0 for matching, 1 for detection, 2+ for nAFC)
%   funccurve: [string] name of the threshold (translation) function
%   parameters: [P+1 unit cell array] where P is the total number of parameters
%       characterizing the funccurve function, and the last parameter is for the
%       psychometric slope.
%       The simplest way to define these parameters is with a cell of 2 units:
%       [v1 v2], which are the bounds of the parameter range
% Output:
%   fast: [struct] FAST structure.
% Example:
%   fast = fastStart(1, 'funcVal', {[1 100], [0.001 1]});
%       for a detection task (nchoice=1)
%       with a constant threshold for all x (funcVal)
%       where that threshold is estimated to be between 1 and 100
%       and the Logistic (default psychometric function) slope is estimated to
%           be between 0.001 and 1
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function fast = fastStart(nchoice, funccurve, parameters)
    fastSettings; %initialize constants.

    funcpsych = psydefault;
    
    fast = fastFull(nchoice, funccurve, funcpsych, parameters);
    
    