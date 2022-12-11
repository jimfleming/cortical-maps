function RetinaRF = make_rf_retina(TotalRetinaCells,RFONCenter,RFOFFCenter,IndRetinaRfspace,onoff_rf_space)
% Make receptive field for each afferent coming to cortex based their
% locationn and polarities in retina

%%
    RetinaRF{TotalRetinaCells} = [];
    filt_rf_retina_on = fspecial('gaussian',60,RFONCenter);  % 19 is the size of the filter it could be width of cortex matrix
    filt_rf_retina_off = fspecial('gaussian',60,RFOFFCenter);
%     filt_rf_retina_on = make_exp_filter(RFONCenter,60);
%     filt_rf_retina_off = make_exp_filter(RFOFFCenter,60);
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


% function filt_ini = make_exp_filter(Xc,size)
%     %Xc = 5; 
%     Xs = round(size/2) ; 
%     %   Difference of Exponential
%     x_interval = 1; 
%     x = 0:x_interval:Xs; 
%     %   Lambda = 0.01 (it means the function has value of 0.01 at x = Xc), 3.0348
%     %   exp( -x^2/2sigma^2) = lambda ==> xc = sigma * sqrt(2*log(1/Lambda))
%     sigma = Xc / 3.034; 
%     %sigma = Xc ; 
%     Yc =  exp( -(x.^2 ) / (2*sigma^2) );
%     DoE = Yc ; 
% 
%     [MeshX,MeshY] = meshgrid(-Xs:Xs,-Xs:Xs);
%     Meshr = round(sqrt( MeshX.^2 + MeshY.^2 )) + 1 ; 
% 
%     filt_ini = zeros( 2*Xs+1 );
%     for i = 1 : Xs + 1 
%         filt_ini( i == Meshr ) = DoE(i);
%     end 
% end 