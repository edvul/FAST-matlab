FAST-matlab
===========

**Functional Adaptive Sequential Testing** toolbox for Matlab

This toolbox is intended for folks conducting multi-trial experiments aimed to estimate thresholds in psychophysics, or wherever else it is necessary. Along the way it does Bayesian estimation of threshold and psychometric functions, thus combining, and extending the functionality of [psignifit](http://psignifit.sourceforge.net/) and [QUEST](http://psych.nyu.edu/pelli/software.html).

**FAST** extends [psignifit](http://psignifit.sourceforge.net/) by jointly estimating the shape of the psychometric function as well as a parametric "threshold function" which describes how the threshold changes with another stimulus variable (e.g., contrast threshold varies as a function of spatial frequency).

Similarly, **FAST** extends classic adaptive testing toolboxes like [QUEST](http://psych.nyu.edu/pelli/software.html) by allowing you to pick stimuli to efficiently, jointly, estimate both the threshold and psychometric functions.  As such, adaptive testing to estimate a single (or multiple independent) thresholds a la _QUEST_ is a special subcase of the functionality of FAST when the threshold function is defined as just a single value.

The details of this method are described in [Vul, Bergsma, & MacLeod (2010)](http://www.evullab.org/pdf/s6.pdf), please cite this if you use FAST in a publication.

## Installation

Copy the code in `source/` to a local directory.  Then in Matlab navigate to that direcoty and run `setup.m` which will add the current directory to the path, allowing you to access FAST functionality.

## Basic Usage

Forthcoming, a more accessible html manual is forthcoming (maybe), in the meantime you can read the [pdf manual](https://github.com/evule/FAST-matlab/raw/master/manual-latex/FAST-manual.pdf).

The basic usage involves a sequence of  
1. initializing the FAST structure by running `fastStart()`  
2. choosing a stimulus value via `fastChooseY()` (or other function)  
3. running a trial with that stimulus value, and updating the FAST structure with the response via `fastUpdate()`  
repeating steps 2 and 3 for a while  
4. extracting threshold and psychometric function estimates via `fastEstimate()`

### Initialization: What you need to know

To use FAST appropriately, you will need to know (and specify) a few properties of the experiment and domain that you are interested in.  This will establish reasonable defaults.  You can also specify stuff with much more precision, but let's leave that for later.

the basic initialization command is: `myfast = fastStart(nchoice, curvefunc, curveparams)`

`nchoice`: **What is the task? (yes/no, matching, nAFC)?**  
This will determine the range, and thus rescaling, of the psychometric function by determining what the level of chance guessing is.  
- `nchoice=0` corresponds to "matching" tasks (tasks where the responses are of the form "too much"/"too little")  
- `nchoice=1` corresponds to "detection" or other "yes"/"no" tasks  
- `nchoice=2+` 1 corresponds to a 2 alternative forced choice task ("this one" or "that one"). `nchoice` can be an arbitrarily large integer for n alternative forced choice tasks with n options.  

`curvefunc`: **What is the threshold function?**  
This will determine how the threshold (y) will vary as a function of some other stimulus variable (x).  The name of the threshold function is provided as a string.  
It is important to consider the domain of stimulus value [linear (-infinity,+infinity) or logarithmic (0, +infinity)], and the appropriateness of the threshold function for that domain.  In general, we find it better to log transform stimulus values (when needed) outside of FAST, and thus stimulus values from the perspective of FAST are linear, so the default psychometric function is a Logistic; but this is not mandatory.  
Different threshold functions included in FAST are:  
- `funcVal`: for estimating a single threshold -- the threshold (y) is just a single value that does not vary with x, and we aim to estimate that value.  
- `funcLine`: for estimating a line -- the threshold (y) varies as a line with x, and we aim to estimate the slope and intercept.  
- `funcPolynomial`: for estimating y as a polynomial function of x -- the number of parameters provided determines the order of the polynomial.  Generally, it is ill-advised to consider polynomial functions of orders much higher than a quadratic.  
- `funcCSF`: for estimating a contrast sensitivity function -- threshold contrast (y) as a function of spatial frequency (x).  There are a few parametric forms of contrast sensitivity functions, and we have a few varities implemented.  
- `funcHyperbolic`: for estimating a simple 1-parameter hyperbolic decay (from 1 to 0 for positive x values)  
- `funcExp`: for estimating an exponential decay/rise function.  Includes three parameters: initial threshold (y) value (when x=0), asymptotic y value (as x >> infinity), and the scaline (e.g., time) constant.  
- `funcTVC`: Threshold vs Contrast function.  (a few varieties are available)  
There are a few other functions built in, but these seemed like the most obvious use cases to start.  It is also reasonably easy to specify your own threshold functions, so long as you make sure they are properly vectorized.  

`curveparams`: **What is a plausible range of parameter values?**  
This is a cell array specifying initial estimates of each parameter of the threshold function, plus the shape of the psychometric function.  
e.g., if we are estimating a simple independent threshold in a logistic psychometric function, we would specify our initial guess about that threshold (say somewhere around 1), and about the slope of the logistic (it should be positive).  We might then specify `curveparams = {[0 2], [0.001 1]}`. The first tuple indicates our initial range for guesses about the threshold, and the second tuple is our initial guess for the logistic slope.  
By default, with a range specified for parameter values, FAST tries to guess whether that parameter should be linear or logarithmic, based on the specified range -- it would assume that the threshold is linear, and the slope is logarithmic here.  FAST then provides a loose normal or log-normal prior on these parameter values given the provided range.  If these defaults are not appropriate, you can specify explicitly whether the parameter should be treated as linear or logarithmic, and what the mean and standard deviation of the normal/log-normal distribution ought to be.

#### Simple example:

If we are using a detection experiment to estimate a single threshold on a linear scale, that ought to be somewhere around 1, we would initialize a FAST structure as:  
`myfast = fastStart(1, 'funcVal', {[0 2], [0.001 1]});`

### 1) Initializaiton

The initialization step amounts to 
