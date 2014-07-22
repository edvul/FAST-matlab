FAST-matlab
===========

**Functional Adaptive Sequential Testing** toolbox for Matlab

This toolbox is intended for folks conducting multi-trial experiments aimed to estimate thresholds in psychophysics, or wherever else it is necessary. Along the way it does Bayesian estimation of threshold and psychometric functions, thus combining, and extending the functionality of [psignifit](http://psignifit.sourceforge.net/) and [QUEST](http://psych.nyu.edu/pelli/software.html).

**FAST** extends [psignifit](http://psignifit.sourceforge.net/) by jointly estimating the shape of the psychometric function as well as a parametric "threshold function" which describes how the threshold changes with another stimulus variable (e.g., contrast threshold varies as a function of spatial frequency).

Similarly, **FAST** extends classic adaptive testing toolboxes like [QUEST](http://psych.nyu.edu/pelli/software.html) by allowing you to pick stimuli to efficiently, jointly, estimate both the threshold and psychometric functions.  As such, adaptive testing to estimate a single (or multiple independent) thresholds a la _QUEST_ is a special subcase of the functionality of FAST when the threshold function is defined as just a single value.

The details of this method are described in [Vul, Bergsma, & MacLeod (2010)](http://www.evullab.org/pdf/s6.pdf), please cite this if you use FAST in a publication.

## Installation

Copy the code to a local directory.  Then in Matlab navigate to that direcoty and run `setup.m` which will add the current directory to the path, allowing you to access FAST functionality.

## Usage

Forthcoming.
