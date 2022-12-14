function [appdata] = make_retina_grid_rf_space_app(appdata);
    rf_space_x = appdata.rf_space_x;
    rf_space_y= appdata.rf_space_y;

    TotalRetinaCells= appdata.TotalRetinaCells;
    RetinaAllX= appdata.RetinaAllX;
    RetinaAllY= appdata.RetinaAllY;
    RetinaAllOD= appdata.RetinaAllOD;
    RetinaAllONOFF= appdata.RetinaAllONOFF;
    RetinotopyRFspace_plot= appdata.RetinotopyRFspace_plot;
    cmap= appdata.cmap;
    rng_trial= appdata.rng_trial;

    % Putting the cells in Matrix  without changing the location (cortical plate)
    % plotting the RGCs on Retinotopic map
    rng(rng_trial)

    TakenSpots = zeros(size(rf_space_x,1),size(rf_space_x,2));
    od_rf_space = zeros(size(rf_space_x));
    onoff_rf_space = zeros(size(rf_space_x));
    IndRetinaRfspace = zeros(size(rf_space_x));

indSelectRand = randperm(TotalRetinaCells);
for ii = 1:TotalRetinaCells
    indTmp = indSelectRand(ii);
    distLocRetCrtx = (rf_space_x(:) - RetinaAllX(indTmp)).^2 + (rf_space_y(:) - RetinaAllY(indTmp)).^2;
    [~,indSelectLoc] = sort(distLocRetCrtx);
    indSelectLoc(TakenSpots(indSelectLoc)==1) = [];
    od_rf_space(indSelectLoc(1)) = RetinaAllOD(indTmp);
    onoff_rf_space(indSelectLoc(1)) = RetinaAllONOFF(indTmp);
    IndRetinaRfspace(indSelectLoc(1)) = indTmp;
    TakenSpots(indSelectLoc(1)) = 1;
end

    appdata.od_rf_space = od_rf_space;
    appdata.onoff_rf_space = onoff_rf_space;
    appdata.IndRetinaRfspace = IndRetinaRfspace;
end
