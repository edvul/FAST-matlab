% function estimate = fastEstimate(fast [,quantileLevels] [,resampleflag] [, quietflag])
% This function estimates parameter quantiles for fast structures.
%
% NOTE: although it is possible to specify parameter sampling densities
% less than 3 in fastInit, fastEstimate does not work with FAST structures
% thusly specified.
%
% Input:
%   FAST structure
% Optional:
%   quantileLevels: by default this is  [.025 .05 .5 .95 .975]
%       that is, the lower and upper bounds of the 95% and 90% CI and the
%       mean.
%   resampleflag: if nonzero (default), start by resampling, to reduce or
%       shift parameter range as appropriate; if zero, don't do that.
%   quietflag: plot (1) or do not plot (0) output
%
% Output:
%   estimate structure:
%   NOTE: on treatment of log parameters
%   Whenever gaussian approximations of parameters are made (the density is
%   estimated as either normal or log-normal) two values are reported:
%   mean and mu.  
%   mu is specified as in the probability distribution;
%   mean is specified linearly.
%   so if a parameter is logarithmic (islog == 1), then 
%   mu == log10(mean). If (islog == 0), then mu == mean
%   sd (standard deviation) is always specified as is natural for the
%   probability distribution: log10() for log-normal distributions.
% 
%   The estimate structure contains estimated quantiles obtained four ways:
%   Parametric Gaussian (estimate.gauss.*): 
%       quadratic interpolation of the marginal distributions for each
%       parameter value.  This yields means and standard deviaions for each
%       parameter, which are combined into a covariance matrix with
%       offdiagonals equal to zero.
%   Empirical Integration and Interpolation (estimate.interp.*):
%       Interpolating the approximated integral of the cross section
%       through the mode.  Unlike all of the above methods, this does not
%       make a gaussian assumption -- in theory these should be most
%       reliable; this needs to be tested.
%   For each of these there will be at least a quantiles cell-array which
%   will contain the parameter values corresponding to the requested or
%   default quantiles.  Also, each will contain parameters of the
%   computation (mu, sigma, and standard deviation for gaussian methods,
%   the parameter values and likelihoods used for integration for
%   interpolation).
% 
% fastEstimate also plots the confidence intervals, marginals and various 
% likelihood surface cross-sections.
%
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function estimate = fastEstimate(fast, quantileLevels, RESAMPLEFLAG, QUIET)

    fastSettings; % load various parameters and settings.

%% set up parameters.
    if nargin < 4
        QUIET = 0;
    end
    if nargin < 3
        RESAMPLEFLAG = 1;
    end    
    if(nargin < 2 || isempty(quantileLevels))
        quantileLevels = defaultQuantiles;
    end
    CIcolors = {'k -'};
    estimate.quantileLevels = quantileLevels;

    % resample?
    GRIDSAMP = cell2mat(fast.params.core.gridsamp);
    if(any(GRIDSAMP == 2))
        error('ERROR: fastEstimate currently does not work with structures with sampling densities < 3');
    end
    if RESAMPLEFLAG
        fast = fastResample(fast, defaultGridSamp, defaultNsigmarange);  
    end


%% Get estimates

    [estimate resamp] = fastCalcEstimates(fast);

%% Compute confidence intervals    
    for i=[1:fast.params.n]
        
        if(fast.params.islog{i})
            estimate.gauss.quantiles{i} = 10.^(norminv(quantileLevels, estimate.gauss.mu(i), estimate.gauss.sd(i)));
%             estimate.marg.quantiles{i} = 10.^(norminv(quantileLevels, estimate.marg.mu(i), estimate.marg.sd(i)));
%             estimate.margXC.quantiles{i} = 10.^(norminv(quantileLevels, estimate.margXC.mu(i), estimate.margXC.sd(i)));
        else
            estimate.gauss.quantiles{i} = (norminv(quantileLevels, estimate.gauss.mu(i), estimate.gauss.sd(i)));
%             estimate.marg.quantiles{i} = (norminv(quantileLevels, estimate.marg.mu(i), estimate.marg.sd(i)));
%             estimate.margXC.quantiles{i} = (norminv(quantileLevels, estimate.margXC.mu(i), estimate.margXC.sd(i)));
        end
    end

%% Compute interpolated confidence intervals...

    maxindex = find(fast.params.core.log10lh(:) == max(fast.params.core.log10lh(:)));
    if(length(maxindex) > 1)
        maxindex = maxindex(ceil(rand()*length(maxindex)));
    end
    [maxindices{1:fast.params.n}] = ind2sub(size(fast.params.core.log10lh),maxindex);
    
    for i = [1:fast.params.n]
        if(range(fast.params.core.pvals{i}) > 0) % is not a singular value (either by weird sampling, or, more likely, by fixing the parameter)

            pts = maxindices;
            pts{i} = [1:length(fast.params.core.pvals{i})];
            llhxc = reshape(squeeze(fast.params.core.log10lh(pts{:})), size(fast.params.core.pvals{i})); % llhxc is log10
            % quantiles by interpolation...cumtrapz for
            % integral of lh; ensure distinctness of values in lhint for interp1
            llhxc = llhxc - max(llhxc);  % max goes to zero.  this is renormalized later and eliminates rounding errors
            lh = 10.^llhxc; % log lh to probability
            lh = lh ./ sum(lh(:)); % normalize.
            if ((lh(1) > max(lh).*10^-2.35) || (lh(end) > max(lh).*10^-2.35)) % ^-4?? dm ^-2.35 is roughly 99.9% confidence interval, no? should be sufficient -ev
                if(~QUIET)
                    fprintf('\nWarning: edges not much less likely than max (param %d): \n\tmaxlh = %f, edges: [%f, %f]\n\t[Especially] interpolated quantiles may be wrong\n\tTry resampling with nsigmas >= 3.3.', i, max(lh), lh(1), lh(end));
                end
            end
            lhint = cumtrapz(fast.params.core.pvals{i},lh);   % rel. fraction of area under lh distn <= x, by approx integration
            %         using trapezoid rule
%             lhint = [cumsum(lhint)];
            lhintN = lhint./max(lhint); % fraction of area under lh distn <= x, by approx integration
            [cumlk, tmpi] = unique(lhintN);
            try
                intQ{i} = interp1(cumlk,fast.params.core.pvals{i}(tmpi),quantileLevels,'cubic');
            catch
                fprintf('\nError interpolating quantiles -- variable dump follows:\n"');
                quantileLevels
                cumlk
                lh
                lhint
                fast.params.core.pvals{i}
                fast.params.core.pvals{i}(tmpi)
            end
        else
            fprintf('\nWarning: singular parameter range for parameter %d\n', i);
            llhxc = repmat(squeeze(fast.params.core.log10lh(maxindices{1})), 1, length(quantileLevels));
            intQ{i} = repmat(fast.params.core.pvals{i}(1), 1, length(quantileLevels));
        end

        % update output structure
        estimate.interp.pvals{i} = fast.params.core.pvals{i}; % the complete vector of (re)sampled parameter i values
        estimate.interp.llhvals{i} = llhxc; % llh10 at those values, with other values set at peak
        estimate.interp.quantiles{i} = intQ{i};
    end

%% plot
if(~QUIET==0)
    QUIET
    figure();
    pvals = fast.params.core.pvals;
    
    % which confidence intervals to plot?
    switch(estimatePLOTCI)
        case {0}
            CIs = estimate.interp.quantiles;
        case {1}
            CIs = estimate.gauss.quantiles;
%         case {2}
%             CIs = estimate.marg.quantiles;
%         case {3}
%             CIs = estimate.margXC.quantiles;
    end
    for i =[1:fast.params.n]
        if(range(fast.params.core.pvals{i}) > 0)
            subplot(fast.params.n, fast.params.n, (i-1)*fast.params.n + i);

            marg = 10.^(fast.params.core.log10lh - max(fast.params.core.log10lh(:))); % convert to probability
            marg = marg./sum(marg(:)); % normalize
            marg = squeeze(sumto1d(marg, i)); % marginalize
            marg = log10(marg); % convert back to log ps
            margll = marg - max(marg(:)); % scale

%             % ## remove later, when cumtrapz works
%             pts = maxindices;
%             pts{i} = [1:length(pvals{i})];
%             llhxc = reshape(squeeze(fast.params.core.log10lh(pts{:})), size(pvals{i})); % llhxc is log10

            plot(pvals{i}, margll, 'b'); % marginal
            hold on;
%             plot(pvals{i}, llhxc-max(llhxc), 'r');  % relative llh10 vs param, XC thru mode

            if(fast.params.islog{i} == 1)
                m = 10.^estimate.gauss.mu(i);
            else
                m = estimate.gauss.mu(i);
            end
            close1 = (pvals{i} - min(CIs{i})).^2;
            close2 = (pvals{i} - max(CIs{i})).^2;
            idx1 = find(close1==min(close1)); % find(min(abs(pvals{i}-min(CIs{i})))
            idx2 = find(close2==min(close2));
            ciylevel = mean([margll(idx1), margll(idx2)])-max(margll);

            axis tight;
            axisdims = axis;
            axisdims(3) = 1.3*(ciylevel-eps); % leave room below ci ylevel
            axis(axisdims);

            plot([m m], [0 axisdims(3)], 'k', 'LineWidth', 2); % Quadratic peak
            xarea = [min(CIs{i}) max(CIs{i})];
            z = fill([xarea fliplr(xarea)], ...
                [0 0 axisdims(3) axisdims(3)], [.3 .3 .3], 'LineStyle', 'none');
            alpha(z, .3);
            plot([min(CIs{i}) max(CIs{i})], [ciylevel ciylevel], 'k:', 'LineWidth', 2);
            axis tight;
            if(i==1)
                legend({'Marginals', 'Quadratic Peak', 'Max CI range'}, 'Location', 'EastOutside');
            end
            if(i<fast.params.n) % current is a curve parameter
                xlabel(sprintf('Curve parameter %d', i));
            else % current is the psych parameters
                xlabel(sprintf('%s slope parameter', func2str(fast.func.psych)));
            end
            ylabel('Log10 likelihood');

            %% plot confidence ellipses...inexact given marginal
            %% approximation
            % ellipses are particularly inexact for planes involving the psy slope
            % parameter, since the contours in that should in principle diverge as
            % slope becomes shallow...but the true contours are checked using
            % contour below, and they seem to agree well enough in practice:
            for j = [i+1:fast.params.n]
                if(range(fast.params.core.pvals{j}) > 0)
                    subplot(fast.params.n, fast.params.n, (j-1)*fast.params.n + i);
                    if(length(fast.params.core.pvals{i}) >= 3 && length(fast.params.core.pvals{j}) >= 3)

                        xp = fast.params.core.pvals{i};
                        yp = fast.params.core.pvals{j};
                        if(fast.params.islog{i})
                            xf = log10(xp);
                        else
                            xf = xp;
                        end
                        if(fast.params.islog{j})
                            yf = log10(yp);
                        else
                            yf = yp;
                        end
                        
                        ll = fast.params.core.log10lh - max(fast.params.core.log10lh(:));
                        Pl = 10.^ll;
                        Pl = Pl ./ sum(Pl(:));
                        for q = [1:fast.params.n]
                            if (q ~= i) & (q ~=j)
                                Pl = sum(Pl,q);
                            end
                        end
                        Pl = squeeze(Pl);
                        
                        LL = log10(Pl);
                        z = LL - max(LL(:));
                        
                        [plotx ploty] = ndgrid(xp, yp);
                        surf(plotx, ploty, double(z));
                        if(defaultColorCutoff<0)
                            caxis([defaultColorCutoff 0]); % color limit cutoff
                        end
                        hold on;
                        contour(plotx,ploty,z, defaultLogLLcontours, 'LineWidth', 2); % outlines redundant with color
                        view(0, 90);
                        % Makes yp or j parameter the vertical axis...
                        v = axis;
                        % prototypical non-elliptical form for yp = psy, xp and yp both
                        % 1:70; if c small, nearly elliptical over many sigma z = ((xp
                        % - mean(mean(xp)))./(yp)).^2  +... c*((yp -
                        % mean(mean(yp)))./(max(max(yp)))).^2


                        [mu covar] = quadEst({xf yf}, z); % 2 and 2x2

                        if((det(covar)>0) && ~(any(diag(covar)<=0)))
                            n=100; % Number of points around ellipse
                            p=pi/n:pi/n:2*pi; % angles around a circle;
                            [eigvec,eigval] = eig(covar); % Compute eigen-stuff
                            xy = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation; eigvals are x,y (co)variances
                            x = xy(:,1);
                            y = xy(:,2);
                            zval = max(x)./sqrt(eigval(1,1));
                            zval2 = max(y)./sqrt(eigval(2,2));
                            % Standard normal deviate along above ellipse contour is
                            % max(x)/sigmax. ? Why not also *exactly* max(y)/sigmay?
                            % Circular standard normal distn has exp(-r^2/2) of its
                            % volume outside radius r, hence at .05 point, r should be
                            % sqrt(5.99) instead of zval; and at .25 point, sqrt(2.77).
                            l = length(quantileLevels);
                            for z = [1:ceil(l/2)]
                                if((quantileLevels(z) + quantileLevels(l+1-z)) == 1) % quantile pair.
                                    p = min(quantileLevels(z), 1-quantileLevels(z))*2; % 2 tail prob
                                    k = sqrt(-2*log(p-.001))/zval; %.001 to eps
                                    if(fast.params.islog{i})
                                        xpl = 10.^((mu(1)+k*x));
                                    else
                                        xpl = (mu(1)+k*x);
                                    end
                                    if(fast.params.islog{j})
                                        ypl = 10.^((mu(2)+k*y));
                                    else
                                        ypl = (mu(2)+k*y);
                                    end
                                    plot(xpl,ypl, CIcolors{mod(z-1, length(CIcolors))+1}, 'LineWidth', 2);
                                    hold on;
                                end
                            end
                            axis tight;
                            if(i<fast.params.n) % current is a curve parameter
                                xlabel(sprintf('Curve param %d', i));
                            else % current is the psych parameters
                                xlabel(sprintf('%s slope param', func2str(fast.func.psych)));
                            end
                            if(j<fast.params.n) % current is a curve parameter
                                ylabel(sprintf('Curve param %d', j));
                            else % current is the psych parameters
                                ylabel(sprintf('%s slope param', func2str(fast.func.psych)));
                            end
                        else
                            fprintf('\nWARNING: Co-Plot %d-%d\nSome sigmas are negative, meaning the log likelihood \nis not properly convex, perhaps we are at a saddle point?\n', i, j);
                        end
                    end
                end
            end
        end
    end
end

