function CxOnOffBal = ONOFF_balance(CXRFON, CXRFOFF, return_abs)
    if nargin < 3
        return_abs = 0;
    end
    tmpCXRFON_norm = max(max(abs(CXRFON)));
    tmpCXRFOFF_norm = max(max(abs(CXRFOFF)));
    CxOnOffBal = (tmpCXRFON_norm - tmpCXRFOFF_norm) / (tmpCXRFON_norm + tmpCXRFOFF_norm);
    if return_abs
        CxOnOffBal = abs(CxOnOffBal);
    end
end