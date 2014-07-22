% function [x,y,(entT)] = fastChooseXYent(fast, nointerp)
%
% Choose (x,y) pair to globally minimize expected
% posterior predictive variance in y* (critical values).  In other words,
% choose a stimulus that is likely to inform the parameter space in a
% manner so as to increase confidence in our predictions about the response
% probability surface. This uses iterative sampling to obtain a good
% estimate quickly. 
% 
%   To use this, you must have initialized fastInit with the xsys parameter
%   specified
% 
% Input:
%   fast: FAST structure
%   Optional
%       nointerp: set to 1 if you don't want to interpolate between XY look-up
%           table values.
% Output:
%   [x y] stimulus pair,
%   Optional
%       eVar: the expected posterior predictive variance
% 
% Example:
%   [x,y] = fastChooseXYpre(MyFAST);
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [varargout] = fastChooseXYpre(fast, evalxs, nointerp)
    if(nargin < 3)
        nointerp = 0;
    end
    if(nargin < 2)
        error('Error: fastChooseYpre must have as input (1) fast (2) a vector of (~10) important x values.');
    end
    
    % ## if we have fixed parameters, quadratic x,y interpolation fails...
    % why? not sure. fix this at some point.
    if any(size(fast.params.core.lattice{1})<3)
        nointerp = 1;
    end
    
    if(fast.params.core.xychoose == 0)
        error('Error: This FAST structure is not configured for fastChooseXY; \nmake sure to run fastInit with an X and Y range specified for this functionality');
    else
        %some setup
        if(length(fast.params.core.xs) < 3)
            xd = 0;
            xm = 1;
        else
            xd = [-1:1];
            xm = 2;
        end
        if(length(fast.params.core.ys) < 3)
            yd = 0;
            ym = 1;
        else
            yd = [-1:1];
            ym = 2;
        end
        
        % find which x and y constraints are in place:
        xidx = [1:length(fast.params.core.xs)];
        yidx = [1:length(fast.params.core.ys)];

        % compute current prior
        bigprior = fast.params.core.log10lh - max(fast.params.core.log10lh(:));
        bigprior = 10.^bigprior;
        bigprior = bigprior ./ sum(bigprior(:));
        bigprior = single(repmat(bigprior+eps, [repmat(1, 1, fast.params.n) length(fast.params.core.xs) length(fast.params.core.ys)]));

        % compute normalizing constants for correct and incorrect.
        normconstC = single(bigprior.*fast.params.core.CPLUT) + eps;
        normconstI = single(bigprior.*(1-fast.params.core.CPLUT)) + eps;
        for i = [1:fast.params.n]
            normconstC = sum(normconstC, i);
            normconstI = sum(normconstI, i);
            idxer{i} = [1:length(fast.params.core.pvals{i})];
            idxconst{i} = 1;
        end

        % compute conditional posterior probabilities
        pvgxc = bigprior.*fast.params.core.CPLUT./repmat(normconstC, [size(fast.params.core.log10lh), 1 1]);
        pvgxi = bigprior.*(1-fast.params.core.CPLUT)./repmat(normconstI, [size(fast.params.core.log10lh), 1 1]);

        % compute for each parameter lattice point, the predicted y for all
        % evalxs x values over which we are maginalizing.
        predys = fastCalcYs(fast, evalxs, -1, fast.params.core.lattice);
        predys = reshape(predys, [numel(fast.params.core.lattice{1}), length(evalxs)]);
        
        BIGYS = repmat(predys, [1 1 numel(fast.params.core.xs) numel(fast.params.core.ys)]);
        
        TC = reshape(pvgxc, [numel(fast.params.core.lattice{1}), 1 numel(fast.params.core.xs) numel(fast.params.core.ys)]);
        BIGTC = repmat(TC, [1, numel(evalxs), 1, 1]);
        varygc = squeeze(mean(sum(BIGYS.^2.*BIGTC,1)./sum(BIGTC,1)-(sum(BIGYS.*BIGTC,1)./sum(BIGTC,1)).^2,2));

        TI = reshape(pvgxi, [numel(fast.params.core.lattice{1}), 1 numel(fast.params.core.xs) numel(fast.params.core.ys)]);
        BIGTI = repmat(TI, [1, numel(evalxs), 1, 1]);
        varygi = squeeze(mean(sum(BIGYS.^2.*BIGTI,1)./sum(BIGTI,1)-(sum(BIGYS.*BIGTI,1)./sum(BIGTI,1)).^2,2));

        normC = squeeze(normconstC);
        normI = squeeze(normconstI);
        
        EVAR = (varygc.*normC + varygi.*normI)./(normC+normI);

        [xi yi] = find(EVAR == min(EVAR(:)));
        xi = xidx(xi(ceil(rand()*length(xi))));
        yi = yidx(yi(ceil(rand()*length(yi))));
        
        % try quadratic interpolation, if it fails, go with lattice min.
        if(nointerp ~= 1)
            xim = max(min(xi, length(fast.params.core.xs)-1), xm);
            yim = max(min(yi, length(fast.params.core.ys)-1), ym);
            
%             [xys{1:2}] = ndgrid(fast.params.core.xs(xd+xim), fast.params.core.ys(yd+yim));
            [muLE sigmaLE] = quadEst2({fast.params.core.xs(xd+xim), fast.params.core.ys(yd+yim)}, -1.*EVAR(xd+xim, yd+yim)); % ND quad interpolation    
            
            if((det(sigmaLE)<=0) || ...
               (muLE(1) < min(fast.params.core.xs(xd+xim))) ||...
               (muLE(1) > max(fast.params.core.xs(xd+xim))) ||...
               (muLE(2) < min(fast.params.core.ys(yd+yim))) ||...
               (muLE(2) > max(fast.params.core.ys(yd+yim))))
                muLE = [fast.params.core.xs(xi) fast.params.core.ys(yi)];
            end
        else
            muLE = [fast.params.core.xs(xi) fast.params.core.ys(yi)];
        end
        
        % output
        if(nargout > 1)
            varargout{1} = double(muLE(1));
            varargout{2} = double(muLE(2));
        else
            error('Error: fastChooseXY must output at least 2 values {x y}, or 3 {x y entTotal}');
        end
        if(nargout == 3)
           varargout{3} = EVAR;
        end
    end
end