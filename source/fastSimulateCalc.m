function fastSimulateCalc(FNAME, plotme, fileout)

% if output file is provided, write to that, else write to screen.
if(nargin<3 || fileout==1)
    fileout = 1;
else
    fileout = fopen(fileout, 'a');
end

% default is to plot
if nargin < 2
    plotme = 1;
end

load(sprintf('%s.mat', FNAME));

NEGLIGIBLELH = .0001;

% if prod(cellfun(@length,TRUTH.P))==1 % sample the parameter space, using final lh, to generate standard errors for y as f(x)
%     allfast = fastResample(allfast); % has the posterior pdf given all data.
%     
%     if(numel(SIM.xs) == 2)
%         evalxs = linspace(min(SIM.xs), max(SIM.xs), 10);
%     else
%         evalxs = SIM.xs;
%     end
%     
%     % compute current prior
%     bigprior = allfast.params.core.log10lh - max(allfast.params.core.log10lh(:));
%     bigprior = 10.^bigprior;
%     bigprior = bigprior ./ sum(bigprior(:));
%     bigprior = single(repmat(bigprior(:)+eps, [numel(bigprior), length(evalxs)]));
%     
%     predys = fastCalcYs(fast, evalxs, -1, fast.params.core.lattice);
%     predys = reshape(predys, [numel(fast.params.core.lattice{1}), length(evalxs)]);
%     
%     meanpredys = sum(predys.*bigprior,1)./sum(bigprior,1);
%     varpredys = var(predys, bigprior, 1);
%     
%     outputs.meanycrit_post_allruns = meanpredys; % the 'marg' mean
%     outputs.stderrycrit_post_allruns =sqrt(varpredys);
% end

fprintf(fileout, '\n\n================================')
fprintf(fileout, '\nRMS and bias in single-run predictions of ycrit, averaged over x and over runs for different estimators')
fprintf(fileout, '\nEst\tRMS err\tbias');
for ee = [1:length(SIM.estimators)]
    outputs.(SIM.estimators{ee}).stdycrit = std(outputs.(SIM.estimators{ee}).ycriterr')'; % sd of predictions, for each x (column)
    outputs.(SIM.estimators{ee}).ycritbias =  mean(outputs.(SIM.estimators{ee}).ycriterr,2); % column, mean error over runs for each x
    outputs.(SIM.estimators{ee}).rmsycriterr =...
        sqrt(sum(outputs.(SIM.estimators{ee}).stdycrit.^2.0...
        + outputs.(SIM.estimators{ee}).ycritbias.^2)./length(SIM.xs)); % scalar
    outputs.(SIM.estimators{ee}).rmsycritbias =  sqrt(mean(outputs.(SIM.estimators{ee}).ycritbias.^2)); % scalar, rms over x of mean error over runs
    fprintf(fileout, '\n%s\t%2.5f\t%2.5f', SIM.estimators{ee}, outputs.(SIM.estimators{ee}).rmsycriterr, outputs.(SIM.estimators{ee}).rmsycritbias);
end

fprintf(fileout, '\n\n================================')
fprintf(fileout, '\nPrecision of estimates of ycrit from all runs combined: comparison to check for overconfidence in posterior pdf')
fprintf(fileout, '\nCol(1), SEM (for combined runs), from (variation in ycrits derived from single-run marg-mean estimates)/sqrt(nruns)')
fprintf(fileout, '\nCol(2), posterior predictive SEM based on overall final pdf from cumulated data: ');
fprintf(fileout, '\n%2.5f\t%2.5f', outputs.marg.stdycrit/sqrt(SIM.nexp-1), 0);%outputs.stderrycrit_post_allruns')

fprintf(fileout, '\n\n================================')
fprintf(fileout, '\nS.D of the distribution of single-run marg-mean final parameter estimates errors');
fprintf(fileout, '\n%2.2f', std(outputs.marg.dev(end,:,:),0,3))
fprintf(fileout, '\nmean S.D. of each run''s final posterior marginal pdf')
fprintf(fileout, '\n%2.2f', mean(outputs.marg.sd(end,:,:),3));  

if(isfield(SIM, 'fxerror') && (SIM.fxerror == 1))
    fast = truefast;
end


%
quantiles = [0.05:0.05:0.95]; % dm fraction falling inside, 2 tailed

for ee = [1:length(SIM.estimators)]
    % mean bias
    aggr.(SIM.estimators{ee}).bias = mean(outputs.(SIM.estimators{ee}).dev, 3);
    aggr.(SIM.estimators{ee}).biasSD = std(outputs.(SIM.estimators{ee}).dev, 0, 3);
    % mean squared error
    aggr.(SIM.estimators{ee}).MSE = mean(outputs.(SIM.estimators{ee}).dev.^2, 3);
    aggr.(SIM.estimators{ee}).MSESD = std(outputs.(SIM.estimators{ee}).dev.^2, 0, 3);
    % L10 squared error
    aggr.(SIM.estimators{ee}).L10DEV = mean(log10(abs(outputs.(SIM.estimators{ee}).dev)+eps), 3);
    aggr.(SIM.estimators{ee}).L10DEVSD = std(log10(abs(outputs.(SIM.estimators{ee}).dev)+eps), 0, 3);
    % posterior SD
    aggr.(SIM.estimators{ee}).L10SD = mean(log10(outputs.(SIM.estimators{ee}).sd+eps), 3);
    aggr.(SIM.estimators{ee}).L10SDSD = std(log10(outputs.(SIM.estimators{ee}).sd+eps), 0, 3);
    % mean Z bias
    aggr.(SIM.estimators{ee}).Z = mean(outputs.(SIM.estimators{ee}).Z, 3);
    aggr.(SIM.estimators{ee}).ZSD = std(outputs.(SIM.estimators{ee}).Z, 0, 3);
    % variance around bias.
    aggr.(SIM.estimators{ee}).var = var(outputs.(SIM.estimators{ee}).dev, 0, 3);
    % confidence interval calibration
    for q = [1:length(quantiles)]
        zthresh = -1.*norminv(0.5-quantiles(q)./2, 0, 1);
        for p = [1:length(datarec{1}.islog)]
            pthresh = sum(sum(abs(outputs.(SIM.estimators{ee}).Z(:,p,:))<=zthresh, 3),1)./numel(outputs.(SIM.estimators{ee}).Z(:,p,:));
            aggr.(SIM.estimators{ee}).QQ(q,p) = pthresh;
        end
    end
end
%%

plotestimators = {'gauss'};%SIM.estimators;

qtry = [0.05 0.1 0.2 0.4 0.6 0.8 0.9 0.95];
qideal = norminv(qtry, 0, 1);
nt = size(aggr.(plotestimators{1}).Z, 1);
np = length(datarec{1}.islog);
if plotme
    colors = {'r-', 'b-', 'g-', 'm-', 'k-'};
    
    % plot the ycrit values
%     figure
%     errorbar(SIM.xs,outputs.meanycrit_post_allruns, outputs.stderrycrit_post_allruns);
%     title('Posterior predictive mean +- SEM of ycrit, based on pdf in param space accumulated over runs')
    
    % plot the across-run average predicted y crit values.
    figure
    plot(SIM.xs, outputs.marg.ycritpredicted, 'r')
    hold on
    plot(SIM.xs, outputs.gauss.ycritpredicted, 'g')
    plot(SIM.xs, T.trueycrit,'b*')
    title('ycrits derived from parameter estimates in each run')
    legend('marg', 'gauss', 'truth')

    figure();
    for p = [1:np]
        subplot(np,1,p);
        plot([1 nt], [0 0], 'k-');
        hold on;
        for ee = [1:length(plotestimators)]
            errorbar([1:nt], aggr.(plotestimators{ee}).bias(:,p), aggr.(SIM.estimators{ee}).biasSD(:,p), colors{mod(ee, length(colors))+1});
            hold on;
        end
        title(sprintf('bias p=%d', p));
    end    
    figure();
    for p = [1:np]
        subplot(np,1,p);
        for ee = [1:length(plotestimators)]
            errorbar([1:nt], aggr.(plotestimators{ee}).MSE(:,p), aggr.(SIM.estimators{ee}).MSESD(:,p), colors{mod(ee, length(colors))+1});
            hold on;
        end
        title(sprintf('MSE p=%d', p));

    end   
    figure();
    for p = [1:np]
        subplot(np,1,p);
        for ee = [1:length(plotestimators)]
            errorbar([1:nt], aggr.(plotestimators{ee}).L10DEV(:,p), aggr.(SIM.estimators{ee}).L10DEVSD(:,p), colors{mod(ee, length(colors))+1});
            hold on;
        end
        title(sprintf('L10DEV p=%d', p));

    end   
    figure();
    for p = [1:np]
        subplot(np,1,p);
        for ee = [1:length(plotestimators)]
            errorbar([1:nt], aggr.(plotestimators{ee}).L10SD(:,p), aggr.(SIM.estimators{ee}).L10SDSD(:,p), colors{mod(ee, length(colors))+1});
            hold on;
        end
        title(sprintf('L10SD p=%d', p));

    end   
    figure();
    for p = [1:np]
        subplot(np,1,p);
        for ee = [1:length(plotestimators)]
            errorbar([1:nt], aggr.(plotestimators{ee}).Z(:,p), aggr.(SIM.estimators{ee}).ZSD(:,p), colors{mod(ee, length(colors))+1});
            hold on;
        end
        title(sprintf('mean Z, Zsdp=%d', p));
    end    
    figure();
        for ee = [1:length(plotestimators)]
        subplot(length(plotestimators),1,ee);
        plot([0 1], [0 1], 'k-', 'Color', [0.3 0.3 0.3], 'LineWidth', 4);
        hold on;
        for p = [1:np]
            plot([quantiles], aggr.(plotestimators{ee}).QQ(:,p), colors{mod(p, length(colors))+1}, 'LineWidth', 2);
            hold on;
        end
        axis square
        title(sprintf('confidence qq p=%d', p));
    end    
    figure();
    for p = [1:np]
        subplot(np,1,p);
        clear qp
        for i = [1:nt]
            qp(i,:,p) = quantile(outputs.gauss.Z(i,p,:), qtry);
        end
        nq = length(qtry);
        for q = [1:ceil(nq/2)]
            plot([1:nt], qp(:,q,p), colors{mod(q, length(colors))+1}, 'LineWidth', 2);
            hold on
            plot([1:nt], qp(:,nq+1-q,p), colors{mod(q, length(colors))+1}, 'LineWidth', 2);
%             plot([1 nt], [qideal(q) qideal(q)], colors{mod(q, length(colors))+1});
%             plot([1 nt], [qideal(nq+1-q) qideal(nq+1-q)], colors{mod(q, length(colors))+1});
        end
        title(sprintf('Z-quant p=%d', p));
        ylim([-2 2])
    end
end

%%
save(sprintf('%s_fin.mat', FNAME), 'datarec', 'SIM', 'fast', 'T', 'outputs' ,'aggr');
