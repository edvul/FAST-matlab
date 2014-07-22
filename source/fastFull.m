% function fast = fastFull(nchoice, funccurve, funcpsych, parameters [, xylims])
% (For an in-depth description of FAST toolbox does, please see the manual)
% fastFull() initializes a FAST structure.
% Input:
%   nchoice: [int] task (0 for matching, 1 for detection, 2+ for nAFC)
%   funccurve: [string] name of the threshold (translation) function
%   funcpsych: [string] name of the psychometric function
%   parameters: [P unit cell array] where P is the total number of params
%       each cell can be a vector of 2, 3, or 4 units.
%       if 2 units [v1 v2] are the bounds of the parameter range
%           this defaults to a uniform prior over the range.
%           e.g., {[0.001 1], ....}
%       if 3 units [S v1 v2]
%           if S=1, assume log-normal prior (mean and sd in log10 scale).
%               In this case we assume that the parameter is referred to
%               by its logarithm. Fast will operate on the logarithm 
%               internally. But the value passed to funccurve can 
%               be the linear value, e.g. 10.^v1, if the funccurve function
%               is written accordingly. Example: funcCSF.m.
%               This convention allows very broad log normal priors, which 
%               are awkward to specify otherwise.
%               v1 is the mean of the decimal log of the parameter value, 
%               v2 is the standard deviation of the log of the value
%           if S=0, assume normal prior (mean and sd in linear scale)
%               v1 is the mean, 
%               v2 is the standard deviation
%           if S=-1, assume uniform prior (bounds on a linear scale)
%               v1 is the lower bound, 
%               v2 is the upper bound,
%               (This is identical to using the 2-unit specification above)
%       if 4 units [S v1 v2 G]
%           same as 3 unit for S, v1, v2, G specifies grid density: 
%           how many points to sample along this parameter (only specify this 
%           if some parameters are more valuable than others, so some should be
%           sampled more finely than others).
%   Optional
%     xsys: [2 unit cell array] Only use this if you plan to use global
%     entropy minimization for stimulus sampling.
%           unit 1 is x, 2 is y.
%           each unit is a 2 unit vector [v1 v2] specifying the upper and
%           lower bounds on the x or y range.
% Output:
%   fast: [struct] FAST structure.
% 
% copyleft Ed Vul & Don MacLeod, 2007
% contact: evul@mit.edu
% version: FAST v0.1

function fast = fastFull(nchoice, funccurve, funcpsych, parameters, xsys)
    fastSettings; %initialize constants.

%% validate input arguments.
    if(nargin < 5)
        xsys = []; % initialize xsys as blank
    end
    if(nargin < 4)
        error(sprintf('Must have at least 4 parameters: \n\tfunccurve\n\tnchoice\n\tfuncpsych\n\tparameters'));
    end
    
    if((~isstr(funccurve)) || (exist(funccurve) == 0))
        error('Funccurve parameter must be a string referring to an existing function');
    end
    if((~isstr(funcpsych)) || (exist(funcpsych) == 0))
        error('Funcpsych parameter must be a string referring to an existing function');
    end

    if((length(parameters{end})>2) && (parameters{end}(1) == 0))
        error('There is no good reason to set the psychometric function slope parameter to be linear (do you really expect it to be negative?)');
    elseif((length(parameters{end})==2) && (min(parameters{end}) <= 0))
        error('Psychometric slope parameter should NOT have a min <=0');
    elseif((length(parameters{end})>2) && (parameters{end}(1) < 0) && (min(parameters{end}(2:3)) <= 0))
        error('Psychometric slope parameter should NOT have a min <=0');
    end
    
%% initialize the basics
    fast.func.curve = str2func(funccurve);
    fast.params.n = numel(parameters);
    fast.data = [];
    fast.func.psych = str2func(funcpsych);
    fast.params.nchoice = nchoice;
    
%% determine sampling requests and resulting density.
    preq = ones(fast.params.n,1);
    pcal = zeros(fast.params.n,1);
    for p = (1:fast.params.n)
        if(numel(parameters{p}) > 3)
            preq(p) = parameters{p}(4);
            if(parameters{p}(4) == 2)
                error('\nERROR: You have indicated that parameter %d have a sampling density of 2.\nParameter densities can be 1 (fixed) or 3+ (for parameters that are not fixed)\nA sampling density of 2 is insufficient for reasonable estimation.', p);
            elseif(parameters{p}(4) == 1 && fastwarnings)
                fprintf('\nParameter %d has a sampling density of 1: it will be fixed.\n', p);
            elseif((parameters{p}(4) < 1) || ~(mod(parameters{p}(4), 1)==0))
                error('\nERROR: Parameter %d has a sampling density of %2.2f: It should be 1 (fixed) or >3', p, parameters{p}(4));
            end
        else
            pcal(p) = 1;
        end
    end
    treq = prod(preq);
    eqsp = min(PARAMGRID_MAX, floor((PARAMGRID_NUMELMAX/treq).^(1/sum(pcal))));
    if(treq > PARAMGRID_NUMELMAX && fastwarnings)
        fprintf('\nWARNING: Total requested parameter lattice size exceeds the preset \nmaximium lattice size (PARAMGRID_NUMELMAX in fastSettings.m) \nthis may cause substantial slowing.');
    elseif(eqsp < PARAMGRID_MIN && fastwarnings)
        fprintf('\nWARNING: Given the total requested parameter lattice, the remaining \nparameters would have sparser sampling than the preset minimum \n(PARAMGRID_MIN fastSettings.m) this may be inefficient');
    end

    preq(pcal==1) = eqsp;

%% make spacing and grid
    for p = (1:fast.params.n)
        if(length(parameters{p}) == 2)
            parameters{p}(2:3) = parameters{p}(1:2);
            parameters{p}(1) = -1;
        end
        fast.params.priorstat{p}.type = parameters{p}(1);
        switch(parameters{p}(1))
            case 1 % parameter is logarithmic, and mean and S.D. of its decimal log are specified
                fast.params.islog{p} = parameters{p}(1);
                fast.params.priorstat{p}.mu = parameters{p}(2);
                fast.params.priorstat{p}.sd = parameters{p}(3);
                fast.params.core.gridsamp{p} = preq(p);
                if(fast.params.core.gridsamp{p} == 1)
                    pvals = fast.params.priorstat{p}.mu;
                else
                    pvals = linspace(...
                       fast.params.priorstat{p}.mu-DEFAULT_RANGE.*fast.params.priorstat{p}.sd,...
                       fast.params.priorstat{p}.mu+DEFAULT_RANGE.*fast.params.priorstat{p}.sd,...
                       fast.params.core.gridsamp{p});
                end
               fast.params.core.pvals{p} = 10.^pvals;
               gridder{p} = -(((pvals - fast.params.priorstat{p}.mu) ...
                               ./(fast.params.priorstat{p}.sd .* sqrt(2))).^2)...
                             ./2.3026;
            case 0 % parameter is linear, with norm mean and S.D. specified
                fast.params.islog{p} = parameters{p}(1);
                fast.params.priorstat{p}.mu = parameters{p}(2);
                fast.params.priorstat{p}.sd = parameters{p}(3);
                fast.params.core.gridsamp{p} = preq(p);
                
                if(fast.params.core.gridsamp{p} == 1)
                    pvals = fast.params.priorstat{p}.mu;
                else
                    pvals = linspace(...
                       fast.params.priorstat{p}.mu-DEFAULT_RANGE.*fast.params.priorstat{p}.sd,...
                       fast.params.priorstat{p}.mu+DEFAULT_RANGE.*fast.params.priorstat{p}.sd,...
                       fast.params.core.gridsamp{p});
                end
               fast.params.core.pvals{p} = pvals;
               gridder{p} = -(((pvals - fast.params.priorstat{p}.mu) ...
                               ./(fast.params.priorstat{p}.sd .* sqrt(2))).^2)...
                             ./2.3026;
            case -1 % parameter is uniform in range... determine if linear or log, no prior
                tmp = [parameters{p}(2) parameters{p}(3)];
                fast.params.core.gridsamp{p} = preq(p);
                if((min(tmp) > 0) && (((abs(diff(log10(tmp)))>ORDERMAGTHRESH)) || (p==fast.params.n))) % if range spans more than 1 orders of magnitude and nothing is below 0 (or it is slope param), do logspace, else linspace
                    fast.params.islog{p} = 1;
                    tmp = log10([min(tmp) max(tmp)]);
                    
                    if(fast.params.core.gridsamp{p} == 1)
                        fast.params.core.pvals{p} = 10.^mean(tmp);
                    else
                        fast.params.core.pvals{p} = 10.^linspace(tmp(1), tmp(2), fast.params.core.gridsamp{p});
                    end
                else
                    fast.params.islog{p} = 0;
                    tmp = [min(tmp) max(tmp)];
                    if(fast.params.core.gridsamp{p} == 1)
                        fast.params.core.pvals{p} = mean(tmp);
                    else
                        fast.params.core.pvals{p} = linspace(tmp(1), tmp(2), fast.params.core.gridsamp{p});
                    end
                end
                fast.params.priorstat{p}.mu = mean(tmp);
                fast.params.priorstat{p}.sd = Inf; % Infinite SD
                gridder{p} = zeros(size(fast.params.core.pvals{p}));
        end
    end
    
    [fast.params.core.lattice{1:fast.params.n}] = ndgrid(fast.params.core.pvals{:});
    
    [lps{1:fast.params.n}] = ndgrid(gridder{:});
    
    fast.params.core.log10lh = zeros(size(lps{1}));
    for p = [1:fast.params.n]
        fast.params.core.log10lh = single(fast.params.core.log10lh + lps{p});
        fast.params.core.lattice{p} = single(fast.params.core.lattice{p});
        mididx = ceil(length(fast.params.core.pvals{p})./2);
%         fast.params.est.initial{p} = fast.params.core.pvals{p}(mididx);
    end
      
    fast = fastUpdate(fast, []); % update all parameter values to current best estimates based on priors.

    if(nargin > 4) % some x an y range was provided.
        XYGRIDSAMP = max(XYGRID_MIN, floor((XYGRID_NUMELMAX./numel(fast.params.core.log10lh))^(1/2))); % how finely to sample look up table?
        if(length(xsys) ~= 2)
            error('Error: XsYs parameter needs to be a 2 unit cell specifying the possible X and Y values.');
        else
            for i = [1:2]
                if(length(xsys{i}) < 1)
                    xsys{i} = 1;
                elseif(length(xsys{i}) == 2)
                    xsys{i} = sort(xsys{i});
                    if((min(xsys{i}) > 0) && ((log10(xsys{i}(2)) - log10(xsys{i}(1)))>ORDERMAGTHRESH))
                        fast.params.core.isxylog{i} = 1;
                        xsys{i} = logspace(log10(xsys{i}(1)), log10(xsys{i}(2)), XYGRIDSAMP);
                    else
                        fast.params.core.isxylog{i} = 0;
                        xsys{i} = linspace(xsys{i}(1), xsys{i}(2), XYGRIDSAMP);
                    end
                elseif(length(xsys{i})>XYGRIDSAMP)
                    if(override)
                        xsys{i} = xsys{i};
                    else
                        idx = round(linspace(1, length(xsys{i}), XYGRIDSAMP));
                        idx = unique(idx);
                        xsys{i} = xsys{i}(idx);
                    end
                else
                    xsys{i} = xsys{i};% ok keep as is.
                end
            end
            fast.params.core.xs = xsys{1};
            fast.params.core.ys = xsys{2};
            
            [bxs bys] = ndgrid(fast.params.core.xs, fast.params.core.ys);
            ycrit = fastCalcYs(fast, reshape(bxs, numel(bxs), 1), -1, fast.params.core.lattice);
            [ps lls] = fastCalcPsLLs(fast, fast.params.core.lattice{fast.params.n}, ycrit, reshape(bys, numel(bys), 1), repmat(1, numel(bys), 1));
            fast.params.core.CPLUT = single(reshape(ps, [size(fast.params.core.lattice{1}), length(fast.params.core.xs), length(fast.params.core.ys)]));
            
            fast.params.core.xychoose = 1;
        end
    else
        fast.params.core.xychoose = 0;
    end
    
    fast.data = [];

    for p = [1:fast.params.n]
        ndp(p) = size(fast.params.core.lattice{1}, p);
    end
    if(fastinittext)
        fprintf('\n====================================');
        fprintf('\nfast Initialized!');
        fprintf('\n\tCurve Function: %s', func2str(fast.func.curve));
        fprintf('\n\tPsychometric Function: %s', func2str(fast.func.psych));
        fprintf('\n\tParameter sampling density: # %d =  %d', [[1:fast.params.n]; ndp]);
        fprintf('\n====================================\n');
    end