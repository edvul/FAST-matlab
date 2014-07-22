% function fastPlot(fast)
%   This plots the fast.data in a FAST structure.
% 
% This provides a plot of:
%   All of the binomial trials as well as the current best function fit.
%   The distribution of p(response) at different values along the X-axis
%   (to determine if we were systematically predicting inappropriate
%   values -- a sign that our assumed functional form is wrong)
%   The psychometric function fit.
%
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function fastPlot(fast)
% load constants
fastSettings;

if(isempty(fast.data))
    error('fastPlot error: no fast.data present in supplied fast structure');
    plotfast.data = 0;
end

%% plot logarithmically or linearly.
if((min(fast.data(:,1))>0) && ...
    ((log10(max(fast.data(:,1))) - log10(min(fast.data(:,1)))) > ORDERMAGTHRESH)) % logspace
    logx = 1;
    xs = 10.^linspace(log10(min(fast.data(:,1))), log10(max(fast.data(:,1))), 100);
else
    logx = 0;
    xs = linspace(min(fast.data(:,1)), max(fast.data(:,1)), 100);
end

% ## make sure this includes all psychometric functions in funcpsych
if(fast.func.psych(0));
    logy = 1;
else
    logy = 0;
end

if(logx && logy)
    fplotxy = str2func('loglog');
    fplotx = str2func('semilogx');
    fploty = str2func('semilogx');
    fploterrx = 'logx';
    fploterry = 'logx';
elseif(logx)
    fplotxy = str2func('semilogx');
    fplotx = str2func('semilogx');
    fploty = str2func('plot');
    fploterrx = 'logx';
    fploterry = 'lin';
elseif(logy)
    fplotxy = str2func('semilogy');
    fplotx = str2func('plot');
    fploty = str2func('semilogx');
    fploterrx = 'lin';
    fploterry = 'logx';
else
    fplotxy = str2func('plot');
    fplotx = str2func('plot');
    fploty = str2func('plot');
    fploterrx = 'lin';
    fploterry = 'lin';
end

%% plot fast.data;
figure();
legendkey = {};
subplot(3,1,1);

% take this out later.
if fast.params.core.xychoose
    myxs = fast.params.core.xs;
    myys = fast.params.core.ys;
    [MX MY] = ndgrid(myxs,myys);
    fplotxy(MX(:), MY(:), 'k*');
    hold on;
end

if(~isempty(fast.data))
    pos = find(fast.data(:,3) >0.5);
    if(~isempty(pos))
        fplotxy(fast.data(pos,1), fast.data(pos,2),'g.');
        legendkey{end+1} = 'R=1';
        hold on;
    end
    neg = find(fast.data(:,3) < 0.5);
    if(~isempty(neg))
        fplotxy(fast.data(neg,1), fast.data(neg,2),'r.');
        legendkey{end+1} = 'R=0';
        hold on;
    end
    neg = find(fast.data(:,3) == 0.5);
    if(~isempty(neg))
        fplotxy(fast.data(neg,1), fast.data(neg,2),'y.');
        legendkey{end+1} = 'R=0.5';
        hold on;
    end
end


colors = {'k -', 'b -', 'c -', 'm -', };
whichp = fastPsyScale(DEFAULT_PLOT_Ps, fast.params.nchoice, 0);
for i=[1:length(whichp)]
    prediction = squeeze(fastChooseYp(fast, xs, whichp(i)));
    fplotxy(xs, prediction, colors{i}, 'LineWidth', 2);
    hold on;
    legendkey{end+1} = sprintf('P = %2.2f', whichp(i));
end
legend(legendkey);
axis tight;
title(sprintf('fast.data and prediction.'));
xlabel('Function X');
ylabel('Function Y');

%%
subplot(3,1,2);

allx = sort(unique(fast.data(:,1)));
for i=[1:length(allx)]
    idx = find(fast.data(:,1) == allx(i));
    px(i) = mean(fast.data(idx,3));
end
fplotx(allx, px, 'b* ');
hold on;
if(length(allx) >= length(fast.data(:,1))./5) % lump fast.data;
    maxx = max(fast.data(:,1));
    minx = min(fast.data(:,1));
    if logx
        grx = maxx/minx;
        tenth = grx.^(1/10);
        indecile = @(x, decile)((((x<(minx.*tenth.^decile+eps)) + (x>=(minx.*tenth.^(decile-1))))));
        centerx = @(i)(10.^mean([log10(minx.*tenth.^i+eps), log10(minx.*tenth.^(i-1))]));
    else
        grx = maxx-minx;
        tenth = grx./10;
        indecile = @(x, decile)((((x<(minx+tenth.*decile+eps)) + (x>=(minx+tenth.*(decile-1))))));
        centerx = @(i)(mean([minx+tenth.*i+eps, (minx+tenth.*(i-1))]));
    end
    

    if(numel(unique(fast.data(:,1))) > 1) % if we do not have just one x value
        for i=[1:10]
            gd = indecile(fast.data(:,1),i);
            idx = find(gd==2);
            x(i) = centerx(i);
            mp(i) = mean(fast.data(idx,3));
            stder(i) = sqrt(1./length(idx));

        end
    else
        x = minx;
        mp = mean(fast.data(:,3));
        stder = sqrt(1./length(fast.data(:,3)));
    end
    ploterr(x, mp, [], stder, 'g*', fploterrx);
end
title('P(response) at Each X');
hold on;
fplotx([min(allx) max(allx)], [.5 .5], 'k -');
ylim([0 1]);
xlim([min(allx)-eps max(allx)+eps]);
xlabel('Function X');
ylabel('P(response)');

%%


yscrit = squeeze(fastCalcYs(fast, fast.data(:,1), -1));

if(logy)
    funcvalues = (fast.data(:,2)./yscrit(:));
    org = logspace(log10(min(funcvalues)), log10(max(funcvalues)), 100);
    rg = org .* yscrit(1);
else
    funcvalues = (fast.data(:,2)-yscrit(:));
    org = linspace(min(funcvalues), max(funcvalues), 100);
    rg = org + yscrit(1);
end
            
pp = fast.func.psych(fast.params.nchoice, ...
                     fast.params.est.('marg').mean{fast.params.n}, ...
                     yscrit(1), rg);

subplot(3,1,3);
fploty(org, pp, 'b -');
hold on;
fulrg = diff([min(funcvalues), max(funcvalues)]);

fploty(funcvalues+(rand(length(funcvalues),1)-.5).*.02.*fulrg, min(max(fast.data(:,3)+(rand(length(funcvalues),1)-.5).*.02, 0), 1),'kx');

tenths = fulrg./10;
if(tenths > 0)
    for i=[1:10]
     gd = (((funcvalues<(min(funcvalues)+tenths.*i+eps)) + (funcvalues>=(min(funcvalues)+tenths.*(i-1)))));
     idx = find(gd==2);
     x(i) = mean([min(funcvalues)+tenths.*i+eps, (min(funcvalues)+tenths.*(i-1))]);
     mp(i) = mean(fast.data(idx,3));
     stder(i) = sqrt(1./length(idx));
    end
else
    x = min(funcvalues);
    mp = mean(fast.data(:,3));
    stder = sqrt(1./length(fast.data(:,3)));
end
ploterr(x, mp, [],stder, 'g*', fploterry);
ylim([0 1]);
xlim([min(funcvalues).*.99 max(funcvalues).*1.01]);
legend({sprintf('%s', func2str(fast.func.psych)), 'Data', 'Lumped data'});
title('Psychometric fit');
if(logy)
    xlabel('Relative Y (Y/ycrit)');
else
    xlabel('Relative Y (Y-ycrit)');
end
ylabel('P(response)');
