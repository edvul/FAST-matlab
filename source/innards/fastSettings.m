% fastSettings
% These are various settings for the FAST toolbox.  Tinker with them if you
% so see fit, but do so with caution.
%
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

% display Matlab warnings?
warning off;

% display FAST warnings?
fastwarnings = 0;

% display initialization text?
fastinittext = 0;


%% Miscellaneous configuration
ORDERMAGTHRESH = 1; 
% to consider something logarithmically spaced it should span this many
% orders of magnitude (as well as having all values > 0)

PARAMGRID_MIN = 7; % minimum (unless otherwise specified) sampling density per parameter
PARAMGRID_MAX = 21; % maximum (if computing, and there are only a few parameters, there is no reason to have a sampling density higher than this)
PARAMGRID_NUMELMAX = 5000; % maximum number of lattice points in parameter lattice

DEFAULT_RANGE = 2; % default sampling range (for fastInit and fastResample) in std. devs: the range will be the mean +/- DEFAULT_RANGE*sigma

% for global entropy minimization, we must define a conditional probability
% lookup table.  This specified the min xy sampling density and max
% sampling density for the entire lookup table (x*y*parameters)
XYGRID_MIN = 9;
XYGRID_NUMELMAX = 150000;

% which estimator to use for roving-range resampling of parameter space?
RESAMPCENTERUSE = 1; % 1:marginal gaussian fit 
RESAMPSDLIM = 7;  

%% Psychometric function settings
% default psychometric function
psydefault = 'psyLogistic';

% default scaled psychometric parameters.
P_nAFC_lapse = .03;         % lapse rate for nAFC expts
P_detect_misslapse = .01;   % detection expts lapse misses
P_detect_FAlapse = .005;    % detection expts lapse hits
P_match_misslapse = .02;    % matching expts lapse misses
P_match_FAlapse = .02;      % matching expts lapse hits
% a "positive" (e.g. 'yes') response.
% These are separate probabilities for inattentive misses and
% for inattentive false alarms: probability of positive response, ppp is
% ppp = FAlapseprob + (1 - misslapseprob-FAlapseprob)*truep
% The inverse function (y from p) is 
% truep = (ppp - FAlapseprob)./(1 - misslapseprob - FAlapseprob)
% For matching experiments, misslapseprob = FAlapseprob seems appropriate.
% lapseprob is less in yes/no than in FC where the wrong button can
% be hit.

%% fastEstimate settings
defaultQuantiles = [.025 .1 .5 .9 .975]; 
%ideally, these should be paired around the mean (.5) to sum to 1.

defaultLogLLcontours = [-0.001, -.35, -.59, -.84, -1.44]; 
% intuitively, these correspond to iso-probability contours, where the
% probability is max_probability*(10^Value), so a Value of -2, would be an
% iso-probability contour at the probability that is 100 times less likely
% than the peak probability.  If the world is gaussian, then the default
% values [-0.001, -.35, -.59, -.84, -1.44] correspond to the MAP estimate,
% the 80% confidence contour, the 90% confidence contour, the 95%
% confidence contour, and the 99% confidence contour.

defaultColorCutoff = -3.3;
% I found it useful to not scale the color map with the graph, but rather
% scale it from 0 to some meaningfully unlikely quantity.  -2.35, in a
% gaussian world is a confidence level of 99.9%.  I figure at that point it
% doesn't make sense to differentiate options, so they will all be colored
% a forboding, unlikely, blue.  If you like more color, perhaps -3.28 will
% be preferable (99.99% confidence interval -- there should rarely be
% anything beyond this, given that fastEstimate resamples to +/- 4 standard
% deviations of the mean).
% Set to 0 if you don't want a fixed color range, but would prefer the
% color-map to change according to the values of the log likelihood surface

% the fastEstimate settings below are more technical, you shouldn't change
% them unless you know what you are doing.
defaultNsigmarange = 3.3; %(+/- sigmas for fastEstimate resampling)
defaultGridSamp = 21; % default sampling density (high, for resolution, since we should be able to afford this offline).

% which confidence interval to plot in fastEstimate?
estimatePLOTCI = 0; %0:interp, 1:gauss, 2:marginal normal, 3:marg through mode normal.

%% fastPlot settings
DEFAULT_PLOT_Ps = [.25 .5 .75]; 
% when plotting the function estimate, we plot some number of probability
% contours -- they are specified here. Note that these are in unscaled
% probabilities, which are then scaled with fastPsyScale into the range
% relevant for the experiment (determined by the nchoice parameter).