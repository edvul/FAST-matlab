% function fast = fastResample(fast, nsamples, nsigmas)
% fastResample resamples the parameter lattice so as to recenter, and, if
% needed, appropriately shrink the spacing.  Thus, with high-dimensional
% parameter spaces, one needn't sacrifice (as much) due to undersampling.
% This re-centers to the marginal means, using marginal SDs as sigma.
% 
% INPUT:
%   fast structure
%   Optional:
%       nsamples, a P-length vector, where P is the number of parameters,
%       or else a scalar to be applied to all parameters. Each entry
%       specifies the new number of sampled values evaluated for each
%       parameter. A value of 0 means accept the existing value from
%       fast.params.core.gridsamp{parameterindex}; a value of 1 means fix this
%       parameter at its best value 
%       nsigmas is likewise either a P-length vector, or a scalar. It is 
%       the number of standard deviations to be spanned by the sampling
%       range in each direction (thus 2 stands for +/- 2sigma). 
%       0 means accept the default from fastSettings().
% OUTPUT:
%   fast structure.
%
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [fast resamp] = fastResample(fast, nsamples, nsigmas)
fastSettings; % load constant settings

USEPRIOR = 1; % by default use prior
QUIET = 1; % dm display details?
RESAMPCENTERUSE = 1;  % dm, for new quadest

% validate input.
if nargin < 3 % nsigmas not set, use default;
    nsigmas = DEFAULT_RANGE.*ones(fast.params.n, 1);
elseif numel(nsigmas) == 1
    nsigmas = nsigmas.*ones(fast.params.n, 1);    
end
nsigmas(nsigmas == 0) = DEFAULT_RANGE;

if nargin < 2 % nsamples not set, use current defaults;
    nsamples = zeros(fast.params.n, 1);
elseif numel(nsamples) == 1
    nsamples = nsamples.*ones(fast.params.n, 1);
end

if(size(fast.data,1) <= 2*(fast.params.n+5) && fastwarnings)
    fprintf('\nWARNING: fastResample might make you veer away from a \nreasonable parameter range if you try to resample \nwith fewer datapoints than (nparam+1)*10\n\n');
end

%% determine sampling requests.
if(sum(nsamples)>0) % reassigning sampling
    preq = ones(fast.params.n,1);
    pcal = zeros(fast.params.n,1);
    for p = (1:fast.params.n)
        if(nsamples(p) > 0)
            preq(p) = nsamples(p);
        else
            if(length(fast.params.core.pvals{p}) ~= 1) % singular value not given
                pcal(p) = 1;  % still zero if singular
            else
                preq(p) = 1;  % still 1 if singular value
            end
        end
    end
    treq = prod(preq);
    eqsp = floor((PARAMGRID_NUMELMAX/treq).^(1/sum(pcal)));
    if(treq > PARAMGRID_NUMELMAX && fastwarnings)
        fprintf('\nWARNING: Total requested parameter lattice size exceeds the preset \nmaximium lattice size (PARAMGRID_NUMELMAX in fastSettings.m) \nthis may cause substantial slowing.');
    elseif(eqsp < PARAMGRID_MIN && fastwarnings)
        fprintf('\nWARNING: Given the total requested parameter lattice, the remaining \nparameters would have sparser sampling than the preset minimum \n(PARAMGRID_MIN fastSettings.m) this may be inefficient');
    end

    preq(pcal==1) = eqsp;
    for i = [1:fast.params.n]
        if((fast.params.core.gridsamp{i} == 1) && (preq(i) > 1))
            if(fastwarnings)
                fprintf('\nWarning: Cannot resample parameter %d from 1 to >1\n', i);
            end
            fast.params.core.gridsamp{i} = 1;
        else
            fast.params.core.gridsamp{i} = preq(i);
        end
    end
end


%% find center
% find marginal mean, std. dev for each parameter.

[est resamp] = fastCalcEstimates(fast);


% assauming we always use marginal quadratic fits.
% 
% switch RESAMPCENTERUSE
%     case {1}
            center.mu = est.gauss.mu;
            center.sd = est.gauss.sd;
%     case {2}
        center.mu = est.marg.mu;
        center.sd = est.marg.sd;
%     case {3}
%         center.mu = est.margXC.mu;
%         center.sd = est.margXC.sd;
% end

%% determine range, sample, compute prior, etc.


for p = [1:fast.params.n] % revise sample spacing:
    
    if(QUIET == 0)
        fprintf('\nParam %d:\t%d\t%0.5g\t%0.5g\n', p, fast.params.islog{p}, est.gauss.mu(p), est.gauss.sd(p)) % dm, else syntax error for 'marginals'
    end
    if(length(fast.params.core.pvals{p})==1)
        % do nothing
        gridder{p} = zeros(1,1);
    else
%         if(fast.params.islog{p})
%             oldpvals(p,:) = reshape(log10(fast.params.core.pvals{p}), 1, numel(fast.params.core.pvals{p}));
%         else
%             oldpvals(p,:) = reshape(fast.params.core.pvals{p}, 1, numel(fast.params.core.pvals{p}));
%         end

        pvals = linspace(...
            center.mu(p)- nsigmas(p).*center.sd(p),...
            center.mu(p)+ nsigmas(p).*center.sd(p),...
            fast.params.core.gridsamp{p});

        if((USEPRIOR == 1) && (fast.params.priorstat{p}.type ~= -1)) % user prior, and prior not uniform
            gridder{p} = -(((pvals - fast.params.priorstat{p}.mu) ...
                ./(fast.params.priorstat{p}.sd .* sqrt(2))).^2)...
                ./2.3026;
        else
            gridder{p} = zeros(size(pvals));
        end

        if(fast.params.islog{p})
            fast.params.core.pvals{p} = 10.^pvals;
        else
            fast.params.core.pvals{p} = pvals;
        end
    end
end

%% build grids, update, etc.
[fast.params.core.lattice{1:fast.params.n}] = ndgrid(fast.params.core.pvals{:});
[lps{1:fast.params.n}] = ndgrid(gridder{:});

fast.params.core.log10lh = zeros(size(lps{1}));
for p = [1:fast.params.n]
    fast.params.core.log10lh = single(fast.params.core.log10lh + lps{p});
    fast.params.core.lattice{p} = single(fast.params.core.lattice{p});
%     mididx = ceil(length(fast.params.core.pvals{p})./2);
end

newdata = fast.data;
fast.data = [];
[fast resamp] = fastUpdate(fast, newdata);

% update the conditional probability look up table to reflect new parameter
% grid.
if(fast.params.core.xychoose == 1)
    if(fast.params.core.isxylog{1})
        xeval = 10.^linspace(log10(min(fast.params.core.xs)), log10(max(fast.params.core.xs)), 100);
    else
        xeval = linspace(min(fast.params.core.xs), max(fast.params.core.xs), 100);
    end
    
    ps = fastPsyScale([0.05 0.95], fast.params.nchoice, 0);
    ymin = min(fastChooseYp(fast, xeval, ps(1)));
    ymax = max(fastChooseYp(fast, xeval, ps(2)));
    
    if(ymin > 0) && ((log10(ymax) - log10(ymin))>ORDERMAGTHRESH)
        fast.params.core.isxyliog{2} = 1;
        fast.params.core.ys = 10.^linspace(log10(ymin), log10(ymax), numel(fast.params.core.ys));
    else
        fast.params.core.isxyliog{2} = 0;
        fast.params.core.ys = linspace(ymin,ymax, numel(fast.params.core.ys));
    end
        
    [bxs bys] = ndgrid(fast.params.core.xs, fast.params.core.ys);
    ycrit = fastCalcYs(fast, reshape(bxs, numel(bxs), 1), -1, fast.params.core.lattice);
    [ps lls] = fastCalcPsLLs(fast, fast.params.core.lattice{fast.params.n}, ycrit, reshape(bys, numel(bys), 1), repmat(1, numel(bys), 1));
    fast.params.core.CPLUT = single(reshape(ps, [size(fast.params.core.lattice{1}), length(fast.params.core.xs), length(fast.params.core.ys)]));
end



