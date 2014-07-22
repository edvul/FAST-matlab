% function [x,y,(entT)] = fastChooseYent(fast, xinput, yinput, iterations)
% Choose y in range yinput for to locally minimize expected posterior 
% entropy given a particular xinput value, this uses iterative sampling to
% obtain a good estimate quickly.
% Input:
%   fast: FAST structure
%   x [scalar]
%   y [2 unit vector] specify the upper and lower bounds of possible ys.
%   Optional
%       iterations: this algorithm starts with the whole y range and
%       iteratively shrinks around the peak.  More iterations means more
%       shrinking, and a finer-grain estimate, defaults to 3, which is fine
%       for most cases.
% Output:
%   [y] stimulus, OR
%   [x y] stimulus coordinates OR
%   [x y entT] with the expected posterior entropy for ys in the final 
%       iteration.
% 
% Example:
%   [x y] = fastChooseY(MyFAST, 10, [.01 1.5], 4);
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [varargout] = fastChooseYent(fast, x, y, iterations)
    fastSettings;
    if(nargin<3)
        error('Error: fastChooseYent must have a scalar x as input and a 2 unit vector (range) for y.');
    elseif(nargin<4)
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
            ys = logspace(log10(min(yrng)), log10(max(yrng)), gridsamp);
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
        end

        % compute conditional posterior probabilities
        pvgxc = bigprior.*LUT./repmat(normconstC, [size(fast.params.core.log10lh), 1 1]);
        pvgxi = bigprior.*(1-LUT)./repmat(normconstI, [size(fast.params.core.log10lh), 1 1]);

        % compute entropies
        lc = log(pvgxc);
        li = log(pvgxi);
        %     toc
        entC = -1.*pvgxc.*lc;
        entI = -1.*pvgxi.*li;
        %     toc
        for j = [1:fast.params.n]
            entC = sum(entC, j);
            entI = sum(entI, j);
        end
        entT = squeeze(entC.*normconstC + entI.*normconstI);
        % combine correct and incorrect entropies
  
        yi = find(entT == min(entT)); % 1D find
        yi = yi(ceil(rand()*length(yi)));
        
        yrng = [ys(max(1, yi-2)), ys(min(gridsamp, yi+2))];
    end
    
    % try quadratic interpolation, if it fails, go with lattice min.
    if(nointerp ~= 1)
        yim = max(min(yi, length(ys)-1), ym);
        ytriple = ys(yd+yim);
        enttriple = entT(yd+yim)';
        c = polyfit(ytriple,enttriple,2);
        cderiv = [2*c(1) c(2)];
        muLE = -c(2)./(2*c(1));  % can only be one; leave sigmas for fastPlot, fastEstimate
        muLE = min(max(muLE,min(ys)),max(ys)); % trim to limits of y range if necessary
        if(muLE ~= min(ys) & muLE ~= max(ys) & c(1) > 0) % minimum entropy point found within range
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
        varargout{3} = entT;
    end
end