% function [est resamp]=fastCalcEstimates(fast)
% Estimate parameters.  Returns Estimate structure:
%
% Mostly used for internals -- no need to interact with it directly...
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [est resamp]=fastCalcEstimates(fast)
    resamp = 0;
%% find MAP
    maxindex = find(fast.params.core.log10lh(:) == max(fast.params.core.log10lh(:)));
    if(length(maxindex) > 1)
        maxindex = maxindex(ceil(rand()*length(maxindex)));
    end
    [maxindices{1:fast.params.n}] = ind2sub(size(fast.params.core.log10lh),maxindex);
    
%% make useful pval struct
    for i = [1:fast.params.n]
        if(fast.params.islog{i})
            pvals{i} = log10(fast.params.core.pvals{i});
        else
            pvals{i} = fast.params.core.pvals{i};
        end
    end
    
    est.islog = fast.params.islog;
    
%% get probabilities from log10lh
    mpt = 10.^(fast.params.core.log10lh - max(fast.params.core.log10lh(:))); % convert log10ll to probability
    mpt = mpt./sum(mpt(:)); % normalize
    
%% Compute marginals 
    [est.marg r1] = calcMarginals(mpt, pvals);
    resamp = min(1, resamp+r1);
    
%% Compute marginal XC through the mode
%      est.margXC = calcMarginalXCs(mpt, maxindices, pvals);
        
%% Compute quadratic approximation of loglh grid around max lattice point.
    [est.gauss r1] = calcGaussian(fast.params.core.log10lh, pvals);
    resamp = min(1, resamp+r1);
    
%% save parameter estimates in their linear designation
    for i=[1:fast.params.n]
        
        if(fast.params.islog{i})
            est.marg.mean{i} = 10.^(est.marg.mu(i));
%             est.margXC.mean{i} = 10.^(est.margXC.mu(i));
            est.gauss.mean{i} = 10.^(est.gauss.mu(i));
        else
            est.marg.mean{i} = (est.marg.mu(i));
%             est.margXC.mean{i} = (est.margXC.mu(i));
            est.gauss.mean{i} = (est.gauss.mu(i));
        end
        est.latticeMax{i} = fast.params.core.pvals{i}(maxindices{i});
    end
end

%%
function [marginals resamp] = calcMarginals(P, pvals)
    resamp = 0;
    for i=[1:length(pvals)]
        Pmarg = reshape(sumto1d(P, i), 1, length(pvals{i}));
        Pmarg = Pmarg ./ sum(Pmarg); % should already be normalized, but why not.
        marginals.mu(i) = sum(Pmarg .* pvals{i}) ./ sum(Pmarg);
        marginals.sd(i) = sqrt(sum(Pmarg .* (pvals{i} - marginals.mu(i)).^2) ./ sum(Pmarg));
        
        if((pvals{i}(1) > (marginals.mu(i) - marginals.sd(i))) || ...
           (pvals{i}(end) < (marginals.mu(i) + marginals.sd(i)))) % if mean is not more than 1 SD away from edge, consider at edge... should expand range or move mean...
            resamp = 1;
        elseif((max(pvals{i})-min(pvals{i}))/marginals.sd(i) > 6) % if lattice is too broad (>6 sds)... 
            resamp = 1;
        end
    end
end
%% is there a signal from calcMarginalsXCs that should signify a need for
%  resampling?  Maybe if the cross section mean is more than 1 sd away from
%  center?  but I don't know if that really means we should resample...
function marginals = calcMarginalXCs(P, maxindices, pvals)
    for i=[1:length(pvals)]
        pts = maxindices;
        pts{i} = [1:length(pvals{i})];
        Pmarg = reshape(P(pts{:}), size(pvals{i}));
        Pmarg = Pmarg ./ sum(Pmarg);
        marginals.mu(i) = sum(Pmarg .* pvals{i}) ./ sum(Pmarg);
        marginals.sd(i) = sqrt(sum(Pmarg .* (pvals{i} - marginals.mu(i)).^2) ./ sum(Pmarg));
    end
end

%%
function [gauss resamp] = calcGaussian(log10lh, pvals)

    [gauss.mu gauss.sigma resamp] = quadEst(pvals, log10lh); % marginal ND quad interpolation

    for i = [1:length(pvals)]
        gauss.sd(i) = gauss.sigma(i,i)^(1/2);
    end
    
    if(any(~isreal(gauss.sd))) 
        resamp = 1;
    end
end