function fastSimulate(FNAME, SIM, TRUTH, SAVEPATH)
% takes three inputs:
% FNAME (name of output file)
% SIM   (structure specifying what the simulated FAST procedure should be)
% TRUTH (structure specifing what the true underling function being
% measured is)
% plotme [0 or 1] indicates whether to plot results
%
% FNAME is a string for the filename.
%
% TRUTH containts...
% .fxy the 'translation/threshold/etc.' function (eg funcExp)
% .fyp the psychometric function (eg psyWeibul)
% .N (nchoice) (description of the task, see FAST help)
% .P a cell array of true parameter values.
%   if an element of P is a single scalar, then that is said to be the true
%   value
%   if an element of p is a 2 unit vector, then a value is chosen uniformly
%   between the min and max of that element to be the true value for the
%   simulation
%   if an element of p is a >2 unit vector, one of those units is chosen to
%   be the true value of that parameter e.g., P{1} = [1 2 3 4 5], then the
%   value of P{1} for a given simulation will be either 1, 2, 3, 4, or 5.
%
% SIM contains the parameters of the simulation:
%   .fast (the parameters for initializing fast
%       .N (nchoice)
%       .fxy (see above)
%       .fyp (see above
%       .P the cell array given to fast to initialize the structure (see
%       fast help for details)
%   .nexp (number of 'experiments'/simulations to run
%   .ntrials (number of trials in an experiment)
%   .ys (plausible y values, if they are to be randomly sampled, or
%   specified for the XY entropy minimization method
%   .xs (plausible x values for the same purpose)
%   .Xplacement (how to choose X)
%     'XYent' jointly pick X and Y to minimize expected posterior entropy
%     'XYpost' jointly pick X and Y to minimize expected variance of the
%               posterior predictive
%     'Xfix' increment through the n possible x values in SIM.xs for
%               situations with a mandatory fixed order of testing
%     'Xrand' choose an x value uniformly between min SIM.xs and max SIM.xs
%   .Yplacement (how to choose Y)
%     'XYent' see above
%     'XYpre' see above
%     'Yiso' method of 1000 staircases (actually n staircases where n =
%               length(SIM.xs)
%     'Ypre' minimize expected variance of posterior predictive
%     'Yp' choose specific P values as given in SIM.ps
%   .ps (for Yplacement = Yp)
%   .fxerror (set to 1 to simulate sampling under the wrong function assumptions.
%       requires that {SIM.truefast}.N .fxy .fyp and .P be set
%   .estimators cell array of names of which 'estimators' provided by FAST
%           are used as estimates to evaluate simulation success.
%
%   I think that may be it....
% Here is a sample of it in use to simulate an exponential decay
%
% FAST SIMULATION TEMPLATE
% FNAME = 'EXP';
%
% % true function
% TRUTH.fxy = 'funcExp';
% TRUTH.fyp = 'psyLogistic';
% TRUTH.P = {[-4 4] [-15 15] [1 50] [.4]};
% TRUTH.N = 0;
%
% SIM.fast.fxy = 'funcExp';
% SIM.fast.fyp = 'psyLogistic';
% SIM.fast.N = TRUTH.N;
% SIM.fast.P = {[0 0 3 1] [0 0 8 7] [1 0.75 .5 7] [1 -1 .8 7]};
%
%
% % these should be paired, so each pair of Xplacement{i} and Yplacement{i}
% % is a single simulation setting
% exps.Xplacement = {'XYpre'};
% exps.Yplacement = {'XYpre'};
%
% SIM.resampling = 1;
% SIM.estimators = {'marg', 'gauss'};
%
% SIM.xs = [1:100];
% SIM.ys = [-20 20];
% SIM.ps = [.2 .5 .8];
%
% SIM.nexp = 1;
% SIM.ntrials = 1000;
%
% for pp = [1:length(exps.Xplacement)]
%     SIM.Xplacement = exps.Xplacement{pp};
%     SIM.Yplacement = exps.Yplacement{pp};
%     fname = sprintf('%s_%s%s', mainname, SIM.Xplacement, SIM.Yplacement);
%
%     fastSimulate(FNAME, SIM, TRUTH);
% end
%
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

% try
    
    fastSettings;
    
    DEBUG = 1; % dm if DEBUG
    
    if (nargin<4)
        SAVEPATH = '.';
    end
    
    for expn = (1:SIM.nexp)
        %%
        T.fxy = str2func(TRUTH.fxy);
        T.fyp = str2func(TRUTH.fyp);
        T.N = TRUTH.N;
        for z = [1:length(TRUTH.P)]
            switch(length(TRUTH.P{z}))
                case {1}
                    T.P{z} = TRUTH.P{z};
                case {2}
                    T.P{z} = abs(diff(TRUTH.P{z})).*rand(1) + min(TRUTH.P{z});
                case {3}
                    if TRUTH.P{z}(1) == 1 % log
                        T.P{z} = 10.^(TRUTH.P{z}(2)+TRUTH.P{z}(3).*randn(1));
                    else
                        T.P{z} = TRUTH.P{z}(2)+TRUTH.P{z}(3).*randn(1);
                    end
                otherwise
                    T.P{z} = TRUTH.P{z}(ceil(rand(1).*length(TRUTH.P{z})));
            end
        end
        if(strcmp(SIM.Yplacement, 'Yiso'))
            fast = fastFull(SIM.fast.N, SIM.fast.fxy, SIM.fast.fyp, SIM.fast.P);
            for xv = [1:length(SIM.xs)]
                if length(SIM.fast.P{end})== 2
                    psychp = [-1 SIM.fast.P{end} 10];
                else
                    psychp = [SIM.fast.P{end}(1:3) 10];
                end
                fasthelp{xv}= fastFull(SIM.fast.N, 'funcVal', SIM.fast.fyp, {[-1 min(SIM.ys) max(SIM.ys) 10] psychp});
            end
        elseif(strcmp(SIM.Xplacement, 'XYent') || strcmp(SIM.Xplacement, 'XYpre'))
            fast = fastFull(SIM.fast.N, SIM.fast.fxy, SIM.fast.fyp, SIM.fast.P, {[min(SIM.xs) max(SIM.xs)]  [min(SIM.ys) max(SIM.ys)]});
        else
            fast = fastFull(SIM.fast.N, SIM.fast.fxy, SIM.fast.fyp, SIM.fast.P);
        end
        
        if(isfield(SIM, 'fxerror') && (SIM.fxerror == 1))
            truefast = fastFull(SIM.truefast.N, SIM.truefast.fxy, SIM.truefast.fyp, SIM.truefast.P);
        end
        
        fprintf('\nExp %d ...', expn);
        tic;
        for i = (1:SIM.ntrials)
            switch(SIM.Xplacement)
                case {'XYent'}
                    [x y] = fastChooseXYent(fast);
                case {'XYpre'}
                    [x y] = fastChooseXYpre(fast, linspace(min(SIM.xs), max(SIM.xs), 10));
                case {'Xfix'}
                    x = SIM.xs(mod(i-1, length(SIM.xs))+1);%rand()*diff(fasting.xs) + min(fasting.xs);
                case {'Xrand'}
                    if(length(SIM.xs) == 2)
                        x = rand()*(max(SIM.xs) - min(SIM.xs)) + min(SIM.xs);
                    else
                        x = SIM.xs(ceil(rand*length(SIM.xs)));
                    end
            end
            switch(SIM.Yplacement)
                case {'XYent'}
                    % do nothing, done above
                case {'XYpre'}
                    % do nothing, done above
                case {'Yent'}
                    y = fastChooseYent(fast, x, SIM.ys);
                case {'Ypre'}
                    y = fastChooseYpre(fast, x, SIM.ys, linspace(min(SIM.xs), max(SIM.xs), 10));
                case {'Yp'}
                    p = SIM.ps(ceil(rand()*(length(SIM.ps))));
                    y = fastChooseYp(fast, x, p);
                case {'Yiso'}
                    jidx = find(SIM.xs == x);
                    p = SIM.ps(ceil(rand()*(length(SIM.ps))));
                    y = fastChooseYp(fasthelp{jidx}, x, p);
            end
            
            ycrit = T.fxy(T.P(1:end-1), x);
            truep = T.fyp(T.N, T.P{end}, ycrit, y);
            r = (rand(1) <= truep);
            %%
            resamp = 0;
            [fast resamp] = fastUpdate(fast, [x y r]);
            if(strcmp(SIM.Yplacement, 'Yiso'))
                [fasthelp{jidx}] = fastUpdate(fasthelp{jidx}, [1 y r]);
            end
            if((SIM.resampling == 1) && ((resamp == 1) || mod(i,30)==0))
                [fast resamp] = fastResample(fast);
            end
            
            
            if(isfield(SIM, 'fxerror') && (SIM.fxerror == 1))
                [truefast resamp] = fastUpdate(truefast, [x y r]);
                if resamp
                    [truefast resamp] = fastResample(truefast);
                end
            end
            
            % estimators...
            if(isfield(SIM, 'fxerror') && (SIM.fxerror == 1))
                for ee = [1:length(SIM.estimators)]
                    datarec{expn}.(SIM.estimators{ee}).mu(i,:) = truefast.params.est.(SIM.estimators{ee}).mu;
                    datarec{expn}.(SIM.estimators{ee}).sd(i,:) = truefast.params.est.(SIM.estimators{ee}).sd;
                end
            else
                for ee = [1:length(SIM.estimators)]
                    datarec{expn}.(SIM.estimators{ee}).mu(i,:) = fast.params.est.(SIM.estimators{ee}).mu;
                    datarec{expn}.(SIM.estimators{ee}).sd(i,:) = fast.params.est.(SIM.estimators{ee}).sd;
                end
            end
            if(DEBUG & length(fast.params.core.pvals) ==2) % mainly to show trial by trial development in fastVal
                contour(fast.params.core.pvals{2}, fast.params.core.pvals{1}, fast.params.core.log10lh, max(max(fast.params.core.log10lh)) + defaultLogLLlevels )
                keyboard % F5 to continue
            end
        end % of trials loop
        if(isfield(SIM, 'fxerror') && (SIM.fxerror == 1))
            datarec{expn}.islog = truefast.params.islog;
        else
            datarec{expn}.islog = fast.params.islog;
        end
        
        datarec{expn}.truth = T;
        
        % compute error, bias, qq confidence.
        
        for p = [1:length(datarec{expn}.islog)]
            islog = datarec{expn}.islog{p};
            if(islog)
                tp = log10(datarec{expn}.truth.P{p});
            else
                tp = datarec{expn}.truth.P{p};
            end
            
            % error, z-score
            for ee = [1:length(SIM.estimators)]
                outputs.(SIM.estimators{ee}).dev(:, p, expn) = datarec{expn}.(SIM.estimators{ee}).mu(:,p)-tp;
                outputs.(SIM.estimators{ee}).Z(:, p, expn) = outputs.(SIM.estimators{ee}).dev(:,p,expn)./datarec{expn}.(SIM.estimators{ee}).sd(:,p);
                outputs.(SIM.estimators{ee}).sd(:, p, expn) = datarec{expn}.(SIM.estimators{ee}).sd(:,p);
            end
        end
        
        yt = eval([SIM.fast.fxy,'(T.P, SIM.xs)']); % true ycrits per x
        T.trueycrit = yt;
        for ee = [1:length(SIM.estimators)]
            yp{ee} = eval([SIM.fast.fxy,'(fast.params.est.',SIM.estimators{ee},'.mean, SIM.xs)']); % predicted ycrits using estimated param vals
            outputs.(SIM.estimators{ee}).ycritpredicted(:,expn) = yp{ee};   % error in prediction of ycrit for this xvalue, run and estimator
            outputs.(SIM.estimators{ee}).ycriterr(:,expn) = yp{ee}-yt;   % column, error in prediction of ycrit for this xvalue, run and estimator
            outputs.(SIM.estimators{ee}).rmsycriterr(expn) = norm(yp{ee}-yt)./sqrt(length(yt));   % scalar, rms error in prediction of ycrit for this run and estimator
        end
%         EV: not sure why we need this...
%             % cumulative fast, meaningful only if all runs are made with same
%             % SIM parameters
%             if prod(cellfun(@length,TRUTH.P))==1
%                 if expn == 1
%                     allfast = fast;
%                 else
%                     allfast.data = [allfast.data; fast.data];
%                 end
%             end
        
        t = toc;
        % change to save more often
        save(sprintf('%s/%s.mat', SAVEPATH, FNAME), 'datarec', 'SIM', 'fast', 'T', 'outputs');
        fprintf('%0.5g', t);
    end % of expn loop
    
    save(sprintf('%s/%s.mat', SAVEPATH, FNAME), 'datarec', 'SIM', 'fast', 'T', 'outputs');
    
% catch
%     errorfile = sprintf('%s_error.mat', FNAME);
%     save(errorfile);
%     fprintf(sprintf('error in fastSimulate: %s\nfull state saved to %s', 'what error?', errorfile));
% end
