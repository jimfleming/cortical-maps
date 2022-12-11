function RetinaRF = make_rf_retina(TotalRetinaCells,RFONCenter,RFOFFCenter,IndRetinaRfspace,onoff_rf_space)
    % Make receptive field for each afferent coming to cortex based their
    % locationn and polarities in retina

    RetinaRF{TotalRetinaCells} = [];
    filt_rf_retina_on = fspecial('gaussian',60,RFONCenter);  % 19 is the size of the filter it could be width of cortex matrix
    filt_rf_retina_off = fspecial('gaussian',60,RFOFFCenter);
    RFONfilt = filt_rf_retina_on ; 
    RFOFFfilt = filt_rf_retina_off ;
    for jj = 1:TotalRetinaCells
        tempGridFilter = IndRetinaRfspace == jj;
        PolarityRFSpace = onoff_rf_space(tempGridFilter); 
        if PolarityRFSpace == 1
            tmpRFON = filter2(RFONfilt,tempGridFilter);
            tmpRFONOFF = (tmpRFON) ./ max(abs(tmpRFON(:)));
        elseif PolarityRFSpace == -1
            tmpRFOFF = filter2(RFOFFfilt,tempGridFilter);
            tmpRFONOFF = (-1*tmpRFOFF) ./ max(abs(tmpRFOFF(:)));
        end
        RetinaRF{jj} = tmpRFONOFF; 
    end
end 
