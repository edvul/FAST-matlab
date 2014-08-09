%  function [mu sigma resamp] = quadEst(params, ll)
%   This function estimates the peaks of the marginal log-likelihood
%   functions for each parameter in turn.
% It uses quadratic interpolation to find the peak of each marginal log
% likelihood, and uses the curvature at the peak to estimate variance. This
% characterization is apt to fail slightly in the tails, so interp is a
% preferable basis for confidence intervals. But sds from this do typically
% agree well with ones directly computed from the marginal pdfs. 

% 1D version, original ND version was quadest2.m
% This 1D version operates on the N marginals independently.
% Although the ND joint peak is not generally the MAP point,
% it may be appropriate for the main intended purpose, which is
% to set center and scale for a resampling of all parameters together.
% Quadest takes, as input, 
% -- PARAMS a cell array of the parameters (the parameter
% lattice) (the points around the Max. LL point)
% -- LL the decimal log likelihood values for each point on that lattice.
%
% As output it provides the N dimensional Mu and Sigma.
% Sigma is a diagonal matrix of the marginal variances.
% Off-diagonal entries of sigma (covariances) are zero in this version.
% Variances are NaN if the maximum marginal is out of the params range.
% In that case, mu is at the appropriate end of the input params range,
% and the resamp flag (third argument, normally zero) is set to suggest
% resampling.

% Singular parameters stuff needed??
% copyleft Ed Vul and Don MacLeod 2007
%   contact: evul@mit.edu

function [mu sigma resamp] = quadEst(pvals, ll)

    % check if any dimensions are singular.
    for i =[1:length(pvals)]
        ds(i) = length(pvals{i});
    end
    if(any(ds<3))
    % rebuild parameter and likelihood grids to remove singular (<3)
    % dimensions
        singulars = 1;
        
        Opvals = pvals;
        clear pvals;

        pref = 1;
        plut = zeros(length(Opvals),1);
        for i = [1:length(Opvals)]
            if(length(Opvals{i}) > 2)
                plut(i) = pref;
                pvals{pref} = Opvals{i}; % using non-singular dimensions
                dims(pref) = length(Opvals{i}); % to reshape ll matrix
                pref = pref + 1;
            else
                ll = mean(ll, i); % for singular (<3) dimensions, take average. -- is this sensible?  This only does anything peculiar if gridsamp for a particular parameter is 2.
            end
        end
        if(numel(dims)==1)
            dims(end+1) = 1;
        end
        ll = reshape(ll, dims);
    else
        singulars = 0;
    end
    nparams = length(pvals);

%% 1D quadratic fit to NEAR-PEAK marginal for each parameter in turn
    sigma= zeros(nparams); % for now, leave covariances zero
    resamp = 0;
    for i =[1:nparams] % marginal relative log10(likelihood):
        marg = 10.^(ll - max(ll(:))); % convert to probability
        marg = marg./sum(marg(:)); % normalize
        marg = squeeze(sumto1d(marg, i)); % marginalize
        marg = log10(marg); % convert back to log ps
        margll(i,:) = marg - max(marg(:)); % scale
        [margmaxll imargmax margc2] = quadmax(margll(i,:)); % constrained within range
        mu(i) = interp1( 1:length(pvals{i}), pvals{i}, imargmax);
        if(imargmax==1 || imargmax == length(pvals{i}))  % or if isNan(c2)
%           if(fastwarnings)
%               fprintf(sprintf('\nWarning: parameter %d:  fit parabola is not concave -- using sample marg mean/sigma', i));
%           end
            sigma(i,i) = (pvals{i}(end)-pvals{i}(1)).^2 /4; % sd set to half the current range
            resamp = 1;
        else
            sigma(i,i) = -0.217/margc2 * (pvals{i}(floor(imargmax+1))-pvals{i}(floor(imargmax))).^2;
        end
    end
%%
    if(singulars)   % if some dimensions were singular....
    % rebuild full dimensionality covariance to fix any singular
    % dimensions.
    % again, add something to skip this if none of the dimensions were
    % singular.
        Nmu = mu;
        Nsigma = sigma;
        clear mu;
        clear sigma;
        for i = [1:length(Opvals)]
            ii = plut(i);
            if(ii == 0)
                mu(i) = mean(Opvals{i}); % this is silly in cases where length(Oparams{i})==2
            else
                mu(i) = Nmu(ii);
            end
            for j = [1:length(Opvals)]
                ij = plut(j);
                if(ii==0 || ij==0)
                    sigma(i, j) = 0;
                else
                    sigma(i,j) = Nsigma(ii,ij);
                end
            end
        end
    end
end
 
function [xmax, ixmax c2] = quadmax(x); 
% find interior local maxima by quadratics fit to triplets of vector x.
% Any ties are resolved by choosing the lowest index, as in Matlab max...

% if max of array is at edge, triplet may suggest max outside range, or
% minimum inside, depending on sign(quadc(1)); other local max may compete.
% Solution based on Damien Garcia, http://www.biomecardio.com: 
xi = x; 
% Fit a quadratic to each triplet of successive parameter values 
% and select the maximum of the within-range extrema 
x1 = x(1:end-2);
x2 = x(2:end-1);
x3 = x(3:end);
xi(2:end-1) = -(x3-x1).^2/8./(x1+x3-2*x2)+x2;

% ... and their corresponding normalized locations
I = zeros(size(xi));
I(2:end-1) = (x1-x3)./(2*x1+2*x3-4*x2);

%% Interpolated maxima...
% extrapolated extrema are ignored
test = I<=-1|I>=1|isnan(xi);
xi(test) = x(test); % replace extrapolations with original cetner 
[xmax,idxmax] = max(xi);
ixmax = I(idxmax) + idxmax;
if (idxmax==1) || (idxmax==length(x))
    c2 = NaN; % anything better??
else
    coeftstmp = polyfit(idxmax-1:idxmax+1, x(idxmax-1:idxmax+1),2);
    c2 = coeftstmp(1);
end
end

 % curvature c2 is multiplier for  (d_idx)^2 in
 % quadratic for x, x drops by 1 if index shifts by sqrt(1/c2)
 % x drops by n if index shifts by sqrt(n/c2)
 % sigma is sqrt(0.5/c2)
