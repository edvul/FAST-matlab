% function [mu sigma] = quadEst(ll, pvals)
%   This function does an N-dimensional quadratic estimation of the
%   log-likelihood function by fitting the log-marginal-likelihood of each
%   parameter with a gaussian (parabola in log space), then combines these
%   into a single ND Gaussian where all off-diagonal covariance values are
%   0 -- thus, there is no parameter correlation in this fit, perhaps
%   makign it somewhat underconfident when parameters are highly
%   correlated, forming a ridge in posterior space.
% 
% this is the marginal version -- the version that fits an ND paraboloid to
% the full log-likelihood surface is quadEst2.m -- although that function
% can account for angled ellipsoidal posteriors better than the independent
% fit can, it can sometimes be unstable when the maximum point is a saddle,
% thus the simpler version is generally preferred -- one should be able to
% sub in quadEst2.m while preserving all other functionality (but I don't
% recommend it).
%
% Provides optional 3rd "resamp" output, if something indicates resampling is
% in order.
% 
% PARAMS a cell array of the parameters (the parameter
% lattice) (the points around the Max. LL point)
% LL the log likelihood values for each point on that lattice.
%
% As output it provides the N dimensional Mu % and Sigma.
% The second output is the covariance matrix of the ND Gaussian
% approximation to the likelihood distribution, implied by adopting
% polymodel for the log10 likelihood....
% Covariances are by definition all zero in this marginal version.
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [mu sigma resamp] = quadEst(pvals, ll)

    printwarnings = 0;

    resamp = 0;
    
    % build large parameter lattice
    nparams = length(pvals);
    
    % compute all the marginals
    likelihood = 10.^(ll-max(ll(:)));
    likelihood = likelihood./sum(likelihood(:));
    for i =[1:nparams]
        margll{i} = log(sumto1d(likelihood,i)+eps); % log base n marginal, nparams triplets as rows)
        pvals{i} = reshape(pvals{i}, length(pvals{i}), 1);
        margll{i} = reshape(margll{i}, length(margll{i}), 1);
        
    end

    sigma= zeros(nparams); % just for now, leave covariances zero
    
    for i =[1:nparams]
        if length(pvals{i})<3
            mu(i) = mean(pvals{i});
            sigma(i,i) = 0;
        else
            quadcoefts{i} = polyfit(pvals{i}, margll{i}, 2);
            
            % catch deviant parabolic fits.
            if any(isnan(quadcoefts{i})) || quadcoefts{i}(1) >= 0
                if(printwarnings)
                    fprintf(sprintf('\nWarning: parameter %d:  fit parabola is not concave -- using sample marg mean/sigma', i));
                end
                marg = exp(margll{i});
                marg = marg./sum(marg(:));
                mu(i) = sum(marg.*pvals{i});
                sigma(i,i) = var(pvals{i}, marg);
            else % all clear
                mu(i) = - 0.5 .* quadcoefts{i}(2) ./ quadcoefts{i}(1);
                sigma(i,i) = -0.5 ./ quadcoefts{i}(1); % diagonal matrix for now
            end
            
            % if the mean is thrown far outside the current range, adopt some kluge heuristics, and throw a warning...
            if ((mu(i)-max(pvals{i})) > (1 * range(pvals{i}))) % far too big
                mu(i) = max(pvals{i});              % restrict to range
                sigma(i,i) = range(pvals{i}).^2;    % set s.d. to be the current range.
                if(printwarnings)
                    fprintf(sprintf('\nWarning: parameter %d is interpolated to be far outside the lattice', i));
                end
                resamp = 1;
            elseif ((mu(i)-min(pvals{i})) < (-1 * range(pvals{i}))) % far too small
                mu(i) = min(pvals{i});              % restrict to range
                sigma(i,i) = range(pvals{i}).^2;    % set s.d. to be the current range.
                if(printwarnings)
                    fprintf(sprintf('\nWarning: parameter %d is interpolated to be far outside the lattice', i));
                end
                resamp = 1;
            end
            
            nsigmarange = (range(pvals{i})./sqrt(sigma(i,i)));
            
            if((nsigmarange > 6))
                resamp = 1;
            end
        end
    end
    
    
    if any(isnan(mu))
        save MYERROR
        error('NANs!')
    end
    
 end
