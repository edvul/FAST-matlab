% function [x,y,(entT)] = fastChooseYpre(fast, x, yrange, evalxs, iterations)
%
% Choose y in range yinput for a given xinput, to locally minimize expected
% posterior predictive variance in y* (critical values).  In other words,
% choose a stimulus that is likely to inform the parameter space in a
% manner so as to increase confidence in our predictions about the response
% probability surface. This uses iterative sampling to obtain a good
% estimate quickly. 
% 
% Input:
%   fast: FAST structure
%   x [scalar]
%   yrange [2 unit vector] the upper and lower bounds of possible ys.
%   evalxs [list] specifying which x values of the posterior predictive to
%       consider for the goal of minimizing the expected posterior 
%       predictive variance. 
%       This -can- be left off, but it is highly discouraged.
%   Optional
%       iterations: this algorithm starts with the whole y range and
%       iteratively shrinks around the peak.  More iterations means more
%       shrinking, and a finer-grain estimate, defaults to 3, which is fine
%       for most cases.
% Output:
%   [y] stimulus, OR
%   [x y] stimulus coordinates OR
%   [x y eVar] with the expected posterior predictive variance for ys in
%       the final iteration.
% 
% Example:
%   [x y] = fastChooseYpre(MyFAST, 10, [.01 1.5], 4);
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [varargout] = fastChooseYpre(fast, x, y, evalxs, iterations)
    fastSettings;
    if(nargin<4)
        error('Error: fastChooseYpre must have as input (1) fast (2) a scalar x (3) a 2 unit vector (range) for y, (4) a vector of (~10) important x values.');
    end
    if(nargin<5)
        iterations = 3;
    end
    if((numel(x)~=1) || (numel(y) ~= 2))
        error('Error: fastChooseY must have a scalar x as input and a 2 unit vector (range) for y.');
    end
    nointerp = 0;
    gridsamp = 9;
    yrng = y;
    
    xd = 0;
    xm = 1;
    yd = [-1:1];
    ym = 2;
    if(min(yrng)<=0)
        islog = 0;
        ys = linspace(min(yrng), max(yrng), gridsamp);
    elseif(abs(diff(log10(yrng)))<ORDERMAGTHRESH)
        islog = 0;
        ys = linspace(min(yrng), max(yrng), gridsamp);
    else
        islog = 1;
        ys = logspace(log10(min(yrng)), log10(max(yrng)), gridsamp);
    end
    
    for i = [1:iterations]
        if(islog)
            ys = 10.^linspace(log10(min(yrng)), log10(max(yrng)), gridsamp);
        else
            ys = linspace(min(yrng), max(yrng), gridsamp);
        end
    
        %compute LUT
        [bxs bys] = ndgrid(x, ys);
        ycrit = fastCalcYs(fast, reshape(bxs, numel(bxs), 1), -1, fast.params.core.lattice);
        [ps lls] = fastCalcPsLLs(fast, fast.params.core.lattice{fast.params.n}, ycrit, reshape(bys, numel(bys), 1), repmat(1, numel(bys), 1));
        LUT = single(reshape(ps, [size(fast.params.core.lattice{1}), 1, length(ys)]));
        
        % compute current prior
        bigprior = fast.params.core.log10lh - max(fast.params.core.log10lh(:));
        bigprior = 10.^bigprior;
        bigprior = bigprior ./ sum(bigprior(:));
        bigprior = single(repmat(bigprior+eps, [repmat(1, 1, fast.params.n) 1 length(ys)]));

        % compute normalizing constants for correct and incorrect.
        normconstC = single(bigprior.*LUT) + eps;
        normconstI = single(bigprior.*(1-LUT)) + eps;
        for j = [1:fast.params.n]
            normconstC = sum(normconstC, j);
            normconstI = sum(normconstI, j);
            idxer{j} = [1:length(fast.params.core.pvals{j})];
            idxconst{j} = 1;
        end
        
        % compute conditional posterior probabilities
        pvgxc = bigprior.*LUT     ./ repmat(normconstC, [size(fast.params.core.log10lh), 1 1]);
        pvgxi = bigprior.*(1-LUT) ./ repmat(normconstI, [size(fast.params.core.log10lh), 1 1]);
        
        % compute for each parameter lattice point, the predicted y
        
        predys = fastCalcYs(fast, reshape(evalxs, numel(evalxs), 1), -1, fast.params.core.lattice);
        predys = reshape(predys, [numel(fast.params.core.lattice{1}), length(evalxs)]);

%         [bxs temp ] = ndgrid(evalxs, [1:10]);
%         predys = fastCalcYs(fast, reshape(bxs, numel(bxs), 1), -1, fast.params.core.lattice);
%         evalys = linspace(min(predys(:)), max(predys(:)), 10);
%         [temp, bys] = ndgrid(evalxs, evalys);
%         
%         [ps lls] = fastCalcPsLLs(fast, fast.params.core.lattice{fast.params.n}, predys, reshape(bys, numel(bys), 1), ones(numel(bys), 1));
%         
%         predps = reshape(ps, [numel(fast.params.core.lattice{1}), length(evalxs)*length(evalys)]);
%         
        for xi = [1:numel(x)]
            for yi = [1:numel(ys)]
                idxer{fast.params.n+1} = xi;
                idxer{fast.params.n+2} = yi;
                idxconst{fast.params.n+1} = xi;
                idxconst{fast.params.n+2} = yi;
                
                pthetac = pvgxc(idxer{:});
                pthetai = pvgxi(idxer{:});
                normC = normconstC(idxconst{:});
                normI = normconstI(idxconst{:});
                
                % computer variance of predys
                varygc = var(predys, pthetac(:)', 1);
                varygi = var(predys, pthetai(:)', 1);
                
                EVAR(xi,yi) = (mean(varygc(:)).*normC + mean(varygi(:)).*normI)./(normC+normI);
            end
        end
        

        % combine correct and incorrect entropies
  
        yi = find(EVAR == min(EVAR)); % 1D find
        yi = yi(ceil(rand()*length(yi)));
        
        yrng = [ys(max(1, yi-2)), ys(min(gridsamp, yi+2))];
    end
    
    % try quadratic interpolation, if it fails, go with lattice min.
    if(nointerp ~= 1)
        yim = max(min(yi, length(ys)-1), ym);
        ytriple = ys(yd+yim);
        enttriple = EVAR(yd+yim);
        c = polyfit(ytriple,enttriple,2);
        cderiv = [2*c(1) c(2)];
        muLE = -c(2)./(2*c(1));  % can only be one; leave sigmas for fastPlot, fastEstimate
        muLE = min(max(muLE,min(ys)),max(ys)); % trim to limits of y range if necessary
        if(muLE ~= min(ys) & muLE ~= max(ys) & c(1) > 0) % minimum posterior predictive variance point found within range
            % do nothing.
        else
            muLE = ys(yi); % can't use interp.
        end
    else
        muLE = ys(yi);
    end
    
    if nargout == 1
        varargout{1} =double( muLE);
    else
        varargout{1} = double(x);
        varargout{2} = double(muLE);  % yLE
    end
    if(nargout == 3)
        varargout{3} = EVAR;
    end
end