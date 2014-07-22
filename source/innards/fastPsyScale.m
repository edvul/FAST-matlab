% function ps=fastPsyCrop(ppp, nchoice, inverse)
% 
% This function rescales psychometric functions.  This should only be used
% by FAST internal functions, and not by the user directly. However, it's
% functionality is described.
%
% The function can operate to scale full ogival (0-1) to some scaled subset
%   inverse == 0
% Or from a scaled function to the full ogival
%   inverse == 1
% 
% The lapse (false positive/negative) rate parameters are loaded from
% fastSettings.
%
% Input:
%   ppp: probability (either from 0 to 1 [unscaled when inverse~=1] or from
%       P(false positives) to 1-P(false negatives) [scaled, when
%       inverse==1]
%   nchoice: 0 = matching; 1 = detection; >1 nAFC
%   inverse: 0: scaling from full ogival; 1: full ogival from scaled P.
% 
% copyleft Ed Vul & Don MacLeod, 2010
% contact: evul@ucsd.edu
% version: 2010-07-16


function ps=fastPsyScale(ppp, nchoice, inverse)

% load default scaling parameters (assumed lapse rates and false P/N rates)
fastSettings;

if(inverse ~= 1) % proceed from full ogival to cropped.
    if(nchoice > 1) % dealing with nAFC (nchoice determines false positive Ps)
        ps = 1./nchoice + (1-P_nAFC_lapse).* ppp.*(1 - 1./nchoice); 
    elseif(nchoice == 1) % dealing with detection (false positives and false negatives are assymetric)
        ps = P_detect_FAlapse + (1 - P_detect_misslapse-P_detect_FAlapse)*ppp; %  correct  for lapses
    else % dealing with matching (false positives and false negatives are symmetric)
        ps = P_match_FAlapse + (1 - P_match_misslapse-P_match_FAlapse)*ppp; %  correct  for lapses
    end
elseif(inverse == 1)  % proceed from cropped, to full ogival function
    if(nchoice > 1) % dealing with nAFC (nchoice determines false positive Ps)
        ps = (ppp-(1./nchoice))./((1-P_nAFC_lapse).*(1-1./nchoice));
        if((min(ps(:)) < eps) || (max(ps(:))>(1-eps)))
            ps = min(max(ps,eps),1-eps); % provide a  warning if this truncation is needed?
            if(fastwarnings)
                fprintf('Warning: Ps provided to psychometric function are asymptotic\nPs should be greater than the false alarm probability, and less than 1-lapse probability \n(%0.5g or %0.5g)', 1./nchoice, 1-P_nAFC_lapse*(1-1/nchoice));
            end
        end
    elseif(nchoice == 1) % dealing with detection (false positives and false negatives are assymetric)
        ps = (ppp - P_detect_FAlapse)./(1 - P_detect_misslapse - P_detect_FAlapse);
             
        if((min(ps(:)) < eps) || (max(ps(:))>(1-eps)))
            ps = min(max(ps,eps),1-eps); % provide a  warning if this truncation is needed?
            if(fastwarnings)
                fprintf('Warning: Ps provided to psychometric function are asymptotic\nPs should be greater than the false alarm probability, and less than 1-lapse probability \n(%0.5g or %0.5g)', 1./nchoice, 1-P_nAFC_lapse*(1-1/nchoice));
            end
        end
    else % dealing with matching (false positives and false negatives are symmetric)
        ps = (ppp - P_match_FAlapse)./(1 - P_match_misslapse - P_match_FAlapse);
        
        if((min(ps(:)) < eps) || (max(ps(:))>(1-eps)))
            ps = min(max(ps,eps),1-eps); % provide a  warning if this truncation is needed?
            if(fastwarnings)
                fprintf('Warning: Ps provided to psychometric function are asymptotic\nPs should be greater than the false alarm probability, and less than 1-lapse probability \n(%0.5g or %0.5g)', 1./nchoice, 1-P_nAFC_lapse*(1-1/nchoice));
            end
        end
    end
end
