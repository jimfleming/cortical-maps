function [od_rf_space,onoff_rf_space,IndRetinaRfspace] = make_retina_grid_rf_space(rf_space_x,rf_space_y,...
    TotalRetinaCells,RetinaAllX,RetinaAllY,RetinaAllOD,RetinaAllONOFF, ...
    RetinotopyRFspace_plot,cmap,rng_trial,show_fig)
% Putting the cells in Matrix  without changing the location (cortical plate)
% plotting the RGCs on Retinotopic map 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    %% ON OFF location in Reina Space
    if show_fig == 1
        figure;
        set(gcf, 'units', 'normalized','outerposition',[0 0 0.9 0.9]);
        subplot(131),imagesc(od_rf_space),axis square; axis off;title('Retina Contra/Ipsi Position','fontsize',16),colormap(gca,'gray'),
        subplot(132),imagesc(onoff_rf_space); colormap(gca,'jet'); axis square; axis off; title('Retina ONOFF Position','fontsize',16), axis off
        subplot(133),imshow(RetinotopyRFspace_plot,cmap); title('Retinotopy','fontsize',16); axis square; axis off

        % Plotting ON OFF RGCs on Retinotopy map
        figure
        set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
        subplot(131),imshow(RetinotopyRFspace_plot,cmap)
        hold on
        NColor = size(rf_space_x,1);
        %   Contra
        MarkerSize = 10;
        CrtxIndRetinaVec = IndRetinaRfspace;
        for p1 = 1 : round(TotalRetinaCells/2)
            [rr1,cc1] = find(CrtxIndRetinaVec == p1);
            Pol1 = onoff_rf_space(rr1,cc1);
            if Pol1 > 0
                plot( cc1 , rr1 ,'r.' ,'MarkerSize', MarkerSize)
            elseif Pol1 < 0
                plot( cc1 , rr1 ,'b.' ,'MarkerSize', MarkerSize)
            end
        end
        title('Contra','fontsize',20)
        xlim([0 NColor]),ylim([0 NColor])

        subplot(132), imshow(RetinotopyRFspace_plot,cmap)
        hold on
        %   Ipsi
        for p1 = round(TotalRetinaCells/2) + 1 : TotalRetinaCells
            [rr1,cc1] = find(CrtxIndRetinaVec == p1);
            Pol1 = onoff_rf_space(rr1,cc1);
            if Pol1 > 0
                plot( cc1 , rr1 ,'r.' ,'MarkerSize', MarkerSize)
            elseif Pol1 < 0
                plot( cc1 , rr1 ,'b.' ,'MarkerSize', MarkerSize)
            end
        end
        title('Ipsi','fontsize',20)
        xlim([0 NColor]),ylim([0 NColor])

        subplot(133), imshow(RetinotopyRFspace_plot,cmap)
        hold on
        %   All Retina on Retinotopy Map
        for p1 = 1 : TotalRetinaCells
            [rr1,cc1] = find(CrtxIndRetinaVec == p1);
            Pol1 = onoff_rf_space(rr1,cc1);
            if Pol1 > 0
                plot( cc1 , rr1 ,'r.' ,'MarkerSize', MarkerSize)
            elseif Pol1 < 0
                plot( cc1 , rr1 ,'b.' ,'MarkerSize', MarkerSize)
            end
        end
        title('Center of Cells On Visuotopy Map','fontsize',20)
        xlim([0 NColor]),ylim([0 NColor])
    end

end
