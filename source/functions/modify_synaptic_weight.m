function [LWON, LWOFF, zero_polarity] = modify_synaptic_weight(rfON, rfOFF, weight_rf, ONOFFCrtxPlt, norm)
    if nargin < 5
        norm = 1;
    end
    
    % Modifying synaptic weights 
    max_resp_on = sum(rfON(:));         % max(allCXrfTempON(:));
    max_resp_off = sum(abs(rfOFF(:)));  % max(abs(allCXrfTempOFF(:)));

    if max_resp_on ~= 0 && max_resp_off ~= 0
        %weight_rf = weight_coef_RF(ii,jj); % .9 to make it lower than the dominant polarity
        if ONOFFCrtxPlt == 1
            LWON = 1;
            if norm
                LWOFF = weight_rf * (max_resp_on/max_resp_off);
            else
                LWOFF = weight_rf;
            end
        elseif ONOFFCrtxPlt == -1
            if norm
                LWON = weight_rf * (max_resp_off/max_resp_on);
            else
                LWON = weight_rf;
            end
            LWOFF = 1;
        end
    else
        LWON = 1;
        LWOFF = 1;
        zero_polarity = 1;  % to measure the number of cortical locations that non-dominant polarity is zero
    end
end
