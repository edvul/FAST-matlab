% function [mu sigma] = quadEst2(params, ll)
%   This function does an N-dimensional quadratic estimation of the
%   log-likelihood function.
% 
% This is based in part on John D'Errico, polyfitn, Mathworks Central file sharing,
% with much sacrifice of generality, some of precision for a little more speed
% This version restores scaling for uniformity in variance during QR
% factoring, though it isn't clear if that makes any difference at all in practice, so
% it could probably be dropped in the interests of speed.
% 
% It takes, as input, 
% PARAMS a cell array of the parameters (the parameter
% lattice) (the points around the Max. LL point)
% LL the log likelihood values for each point on that lattice.
%
% Provides optional 3rd "resamp" output, if something indicates resampling is
% in order.
% 
% As output it provides the N dimensional Mu % and Sigma.
% The third output d2quadvaldpi gives the partial second derivatives as column vectors
% in an nparams by npts matrix, where the ith row is the derivative wrt
% the ith parameter value. This determines the covariance matrix...
% The 4th argument is the covariance matrix of the ND Gaussian
% approximation to the likelihood distribution, implied by adopting
% polymodel for the log10 likelihood....
% 
% copyleft Ed Vul and Don MacLeod 2007
%   contact: evul@mit.edu

function [mu sigma resamp] = quadEst2(pvals, ll)

    resamp = 0;
    
    nparams = length(pvals);
    
    % Moved this machinery into quadEst2 from calcEstimates to maintain
    % generality with quadEst while quadEst calculates full marginals.
    maxindex = find(ll(:) == max(ll(:)));
    if(length(maxindex) > 1)
        maxindex = maxindex(ceil(rand()*length(maxindex)));
    end
    [maxindices{1:nparams}] = ind2sub(size(ll),maxindex);
    
    pts = maxindices;
    for i =[1:nparams]
        if(length(pvals{i}) < 3)
            pv{i} = pvals{i}(maxindices{i});
        else
            pts{i} = min(max(pts{i},2), length(pvals{i})-1);
            pts{i} = max(pts{i}, 2);
            if(maxindices{i} ~= pts{i}) % if MAP was at edge and we moved, likely should resample.
                resamp = 1;
            end
            pts{i} = pts{i} + [-1:1];
            pv{i} = pvals{i}(pts{i});
        end
    end
    ll = ll(pts{:}) .* 2.3026; % convert to natural log while we are at it.
    

    % check if any dimensions are singular.
    for i =[1:length(pvals)]
        ds(i) = length(pvals{i});
    end
    if(any(ds<3))
    % rebuild parameter and likelihood grids to remove singular (<3)
    % dimensions
        singulars = 1;
        
        Opvals = pvals;
        Oll = ll;
        clear pvals;

        pref = 1;
        plut = zeros(length(Opvals),1);
        for i = [1:length(Opvals)]
            if(length(Opvals{i}) > 2)
                plut(i) = pref;
                pvals{pref} = Opvals{i}; % using non-singular dimensions
                dims(pref) = length(Opvals{i}); % to reshape ll matrix
                pref = pref + 1;
            else
                ll = mean(ll, i); % for singular (<3) dimensions, take average. -- is this sensible?  This only does anything peculiar if gridsamp for a particular parameter is 2.
            end
        end
        if(numel(dims)==1)
            dims(end+1) = 1;
        end
        ll = reshape(ll, dims);
    else
        singulars = 0;
    end
    nparams = length(pvals);
    if(nparams == 1)
        params{1} = pvals{:};
    else
        [params{1:nparams}] = ndgrid(pvals{:});
    end
    depvar = ll(:);
    indepvar = [];
    for i =[1:nparams]
        indepvar(:,i) = params{i}(:);
    end
    polymodel.termexponents = polynbasis(2,nparams);
    nt = size(polymodel.termexponents,1); % number of terms in quadratic model
    ndata = length(depvar);

    stdind = sqrt(diag(cov(indepvar)));
    % scaled variables
    indepvar_s = indepvar*diag(1./stdind);
    % build the design matrix
    M = ones(ndata,nt);
    scalefact = ones(1,nt);
    for i = 1:nt
      for j = 1:nparams
        M(:,i) = M(:,i).*indepvar_s(:,j).^polymodel.termexponents(i,j);
        scalefact(i) = scalefact(i)/(stdind(j)^polymodel.termexponents(i,j));
      end
    end

    %     M = ones(ndata,nt);
    %     for i = 1:nt
    %         for j = 1:nparams
    %             M(:,i) = M(:,i).*indepvar(:,j).^polymodel.termexponents(i,j);
    %         end
    %     end
    
    % (D'Errico): estimate the model using QR. do it this way to provide a
    % covariance matrix when all done. Use a pivoted QR for
    % maximum stability.
    [Q,R,E] = qr(M,0);
    rhs = Q'*depvar;
    polymodel.Coefficients(E) = R\rhs;
    yhat = M*polymodel.Coefficients';
    
    polymodel.Coefficients=polymodel.Coefficients.*scalefact; % M and scaling are history now; omit if no scaling above

    % construct coefficient matrix for partial derivatives of poly
    % row i has coefficients to express d(poly)/dpval(i), as wtd sum of pvals
    % Thus, nparams by nparams matrix with one coefficient for each pair of
    % parameters

    coefts = polymodel.Coefficients'; % nt by 1
    [i2,j2] = find(polymodel.termexponents == 2);
    pderivterms=diag(2*coefts(i2)); % diagonal terms from squared components

    nsi = setdiff(1:nt,i2);   % index to remaining (non-squared) components in termexponents, coefts

     crossprods = find(sum(polymodel.termexponents(nsi,:)') == 2);
    cpi = nsi(crossprods); % cpi indexes the crossproduct rows in termexponents and coefts

    % termexponents(cpi,:) is npar(npar-1)/2 by npar matrix of crossp exponents

    for ii = 1:length(cpi)
        [rowparam, colparam] = find(polymodel.termexponents(cpi(ii),:) > 0);
        pderivterms(colparam(1), colparam(2)) = coefts(cpi(ii));
        pderivterms(colparam(2), colparam(1)) = coefts(cpi(ii));  % lower triangle
    end
    lineari = setdiff(nsi,cpi); % index to npar+1 linear terms, inc final constant
    lineari = lineari(1:end-1); % assume the constant term is indeed final...
    rhsvals = -coefts(lineari); % rhsvals(i) is constant in pderivterms for param(i), from linear term
    mu = pderivterms\rhsvals;  % column vector of parameter values at maximum; 
    % first partial derivatives at mu, pderivterms*mu-rhsvals, should be nearly zero:
    % pderivterms*mu-rhsvals  % debug check  
    % Conveniently, the matrix of second-order partial derivatives is just pderivterms/2; inversion gives cov matrix:
    sigma = -(inv(pderivterms));   % ## /2? dm, previously not /2; Did this ever matter?
    % second-order partial derivatives=pderiv/2
    % the comments above don't make sense to me -- I think these are from a
    % few revisions back?
    

    
    if(singulars)   % if some dimensions were singular....
    % rebuild full dimensionality covariance to fix any singular
    % dimensions.
    % again, add something to skip this if none of the dimensions were
    % singular.
        Nmu = mu;
        Nsigma = sigma;
        clear mu;
        clear sigma;
        for i = [1:length(Opvals)]
            ii = plut(i);
            if(ii == 0)
                mu(i) = mean(Opvals{i}); % this is silly in cases where length(Oparams{i})==2
            else
                mu(i) = Nmu(ii);
            end
            for j = [1:length(Opvals)]
                ij = plut(j);
                if(ii==0 || ij==0)
                    sigma(i, j) = 0;
                else
                    sigma(i,j) = Nsigma(ii,ij);
                end
            end
        end
    end
            
    
 end


% function based on John D'Errico, Matlab  central file sharing:
% can be run ahead of data input
function polynterms = polynbasis(order, p)
% build the exponent array recursively
% arguments: (input)
%  p: number of parameters (indep vars)
%
% arguments: (output)
%  polynterms - exponent array for the model, nterms by p

    if p == 0
        % terminal case
        polynterms = [];
    elseif (order == 0)
        % terminal case
        polynterms = zeros(1,p);
    elseif (p==1)
        % terminal case
        polynterms = (order:-1:0)';
    else
        % general recursive case
        polynterms = zeros(0,p);
        for k = order:-1:0
            t = polynbasis(order-k,p-1);
            nt = size(t,1);
            polynterms = [polynterms;[repmat(k,nt,1),t]];
        end
    end

end