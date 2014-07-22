% function [x,y,(entT)] = fastChooseXYent(fast, nointerp)
% Choose (x,y) pair to globally minimize expected posterior entropy.
%   To use this, you must have initialized fastInit with the xsys parameter
%   specified, to build the conditional probability look-up table at init
%   initialization (required to make these calculations reasonably quick)
% 
% Input:
%   fast: FAST structure
%   Optional
%       nointerp: set to 1 if you don't want to interpolate between XY look-up
%           table values.
% Output:
%   [x y] stimulus pair,
%   Optional
%       entT: the expected entropy for each x,y grid point
% 
% Example:
%   ys = fastChooseXYent(MyFAST, [1:10], 0.75);
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16

function [varargout] = fastChooseXYent(fast, nointerp)
    if(nargin < 2)
        nointerp = 0;
    end
    if(nargin < 1)
        error('Error: FAST structure needed as input for fastChooseXY');
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
        % find which x and y constraints are in place:
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
        
        xidx = [1:length(fast.params.core.xs)];
        yidx = [1:length(fast.params.core.ys)];

        % compute current posterior (prior for the next trial)
        bigprior = fast.params.core.log10lh - max(fast.params.core.log10lh(:));
        bigprior = 10.^bigprior;
        bigprior = bigprior ./ sum(bigprior(:));
        bigprior = single(repmat(bigprior+eps, [repmat(1, 1, fast.params.n) length(fast.params.core.xs) length(fast.params.core.ys)]));

        % compute normalizing constants for correct and incorrect.
        normconstC = single(bigprior.*fast.params.core.CPLUT) + eps;
        normconstI = single(bigprior.*(1-fast.params.core.CPLUT)) + eps;
        
        % marginalize over all parameters
        for i = [1:fast.params.n]
            normconstC = sum(normconstC, i);
            normconstI = sum(normconstI, i);
        end

        % compute conditional posterior probabilities
        pvgxc = bigprior.*fast.params.core.CPLUT./repmat(normconstC, [size(fast.params.core.log10lh), 1 1]);
        pvgxi = bigprior.*(1-fast.params.core.CPLUT)./repmat(normconstI, [size(fast.params.core.log10lh), 1 1]);

        % compute entropies
        lc = log(pvgxc);
        li = log(pvgxi);
        entC = -1.*pvgxc.*lc;
        entI = -1.*pvgxi.*li;
        % entropy is sum over parameters.
        for i = [1:fast.params.n]
            entC = sum(entC, i);
            entI = sum(entI, i);
        end
        % combine correct and incorrect entropies
        entT = entC.*normconstC + entI.*normconstI;
        entT = reshape(entT, length(fast.params.core.xs), length(fast.params.core.ys));

        entfind = entT(xidx, yidx);
        [xi yi] = find(entfind == min(entfind(:)));
        xi = xidx(xi(ceil(rand()*length(xi))));
        yi = yidx(yi(ceil(rand()*length(yi))));
        
        % try quadratic interpolation, if it fails, go with lattice min.
        if(nointerp ~= 1)
            xim = max(min(xi, length(fast.params.core.xs)-1), xm);
            yim = max(min(yi, length(fast.params.core.ys)-1), ym);
            
%             [xys{1:2}] = ndgrid(fast.params.core.xs(xd+xim), fast.params.core.ys(yd+yim));
            [muLE sigmaLE] = quadEst({fast.params.core.xs(xd+xim), fast.params.core.ys(yd+yim)}, -1.*entT(xd+xim, yd+yim)); % ND quad interpolation    
            
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
%             [varargout{1} varargout{2}]
        else
            error('Error: fastChooseXY must output at least 2 values {x y}, or 3 {x y entTotal}');
        end
        if(nargout == 3)
           varargout{3} = entT;
        end
    end
end