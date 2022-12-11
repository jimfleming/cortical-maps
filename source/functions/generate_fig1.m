function generate_fig1(length_crtx_present,point_x,point_y,ODCrtxPlt_smooth,ONOFF_smoothed,rf_space_x,...
    Retinotopy3mIndex,RetinaAllX,RetinaAllY,RetinaAllOD,RetinaAllONOFF,Retinotopy3mIndPlot,RetinotopyRFspace_plot,RFONCenter,show_fig)
    % parameters to generate figure 1
    % length_crtx_present        % the number of afferents shown in  figure1
    % point_x                    % starting col for the selected area  
    % point_y                    % starting row for the selected area  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if show_fig ==1 
        rect_x1 = point_x - .5 ;
        rect_x2 = rect_x1 + length_crtx_present;
        rect_y1 = point_y - .5 ;
        rect_y2 = rect_y1 + length_crtx_present;
        row_range_fig1 = point_y : (point_y+length_crtx_present-1);
        col_range_fig1 = point_x : (point_x+length_crtx_present-1);

        [ret_map_fig1,cmap_fig1]= show_RGC_colormap(rf_space_x,rect_x1,rect_x2,...
                                            rect_y1,rect_y2,Retinotopy3mIndex,RetinaAllX,RetinaAllY,RetinaAllOD,RetinaAllONOFF,1);

        show_onoff_od_retinotopy(ODCrtxPlt_smooth,ONOFF_smoothed,row_range_fig1,col_range_fig1,...
            length_crtx_present,Retinotopy3mIndPlot,RetinotopyRFspace_plot,ret_map_fig1,cmap_fig1) 
        
        show_retinotopy_with_receptive_field(row_range_fig1,col_range_fig1,length_crtx_present,Retinotopy3mIndPlot,RetinotopyRFspace_plot,ret_map_fig1,cmap_fig1,RFONCenter)
    end 
    
end 

%% functions 
function [ret_map_fig1,cmap_fig1]= show_RGC_colormap(rf_space_x,rect_x1,rect_x2,rect_y1,rect_y2,...
                Retinotopy3mIndex,RetinaAllX,RetinaAllY,RetinaAllOD,RetinaAllONOFF,show_fig)
    % show the location of afferents on retinotpic colormap for the selected area of cortex (figure1 first row)   
    
    index_aff_select = Retinotopy3mIndex(floor(rect_y1:rect_y2),floor(rect_x1:rect_x2));
    ret_x_aff_select = zeros(1,length(index_aff_select(:)));
    ret_y_aff_select = zeros(1,length(index_aff_select(:)));
    retina_od = zeros(1,length(index_aff_select(:)));
    retina_onoff = zeros(1,length(index_aff_select(:)));
    for ww = 1 : length(index_aff_select(:))
        ind = index_aff_select(ww); 
        ret_x_aff_select(ww) = RetinaAllX(ind);
        ret_y_aff_select(ww) = RetinaAllY(ind);
        retina_od(ww) = RetinaAllOD(ind);
        retina_onoff(ww) = RetinaAllONOFF(ind);
    end 
    
    row_min = round(min(ret_y_aff_select(:)) - 1); 
    row_max = round(max(ret_y_aff_select(:)) + 1); 
    col_min = round(min(ret_x_aff_select(:)) - 1); 
    col_max = round(max(ret_x_aff_select(:)) + 1); 
    row_length = row_max - row_min; 
    col_length = col_max - col_min;
    ret_box_length = max(row_length,col_length)+6; % the shift is because to avoid touching to the borders 

    NColor = ret_box_length;
    Retinotopy_map_plot_fig1 = zeros(NColor,NColor);
    Retinotopy_map_plot_fig1(:) = 1 : NColor * NColor ;
    [RedGrid,~] = meshgrid( linspace(0,255,NColor) , linspace(0,0,NColor) );
    [~,GreenGrid] = meshgrid( linspace(230,230,NColor) , linspace(230,51,NColor) );
    [BlueGrid,~] = meshgrid( linspace(0,255,NColor) , linspace(0,0,NColor) );
    cmap_fig1 = [RedGrid(:) GreenGrid(:) BlueGrid(:)]/255;

%     figure, %imshow(Retinotopy_map_plot_fig1,cmap_fig1)
%     imagesc(Retinotopy_map_plot_fig1)
%     colormap(cmap_fig1)
%     axis square

    ret_map_fig1 = nan(size(rf_space_x)); 
    row_range_box = row_min-3:(row_min+ret_box_length-4);
    col_range_box = col_min-3:(col_min+ret_box_length-4);
    ret_map_fig1(row_range_box,col_range_box) = Retinotopy_map_plot_fig1 ;  %  shifted to no touch the borders 
    %% showing all afferents within the selected box on colormap
    
    ind_od_onoff = zeros(length(RetinaAllX(:)),1); 
    for ww2 = 1 : length(RetinaAllX(:))
        ind = ww2;
        if (RetinaAllX(ind)>col_min && RetinaAllX(ind)<col_max && ...
                RetinaAllY(ind)>row_min && RetinaAllY(ind)<row_max)
            if RetinaAllOD(ind) == 1 
                if RetinaAllONOFF(ind) == 1
                    ind_od_onoff(ww2) = 1; 
                elseif RetinaAllONOFF(ind) == -1
                    ind_od_onoff(ww2) = 2; 
                end
            elseif RetinaAllOD(ind) == -1 
                if RetinaAllONOFF(ind) == 1
                    ind_od_onoff(ww2) = 3; 
                elseif RetinaAllONOFF(ind) == -1
                    ind_od_onoff(ww2) = 4; 
                end
            end 
        end 
    end

    ret_map_fig2 = ret_map_fig1; 
    %ret_map_fig2(130,44:55) = 37 ; % for scale bar 
    if show_fig == 1 
        MarkerSize = 60;
        figure, %imshow(Retinotopy_map_plot_fig1,cmap_fig1)
        subplot(121)
        imagesc(ret_map_fig2)
        colormap(cmap_fig1)
        axis square
        hold on
        plot( RetinaAllX(ind_od_onoff==1) , RetinaAllY(ind_od_onoff==1) ,'r.' ,'MarkerSize', MarkerSize)
        plot( RetinaAllX(ind_od_onoff==2), RetinaAllY(ind_od_onoff==2),'b.' ,'MarkerSize', MarkerSize)
        title('Contra','fontsize',20)
        xlim([col_range_box(1) col_range_box(end)]),ylim([row_range_box(1) row_range_box(end)])
        axis on 
        
        subplot(122)
        imagesc(ret_map_fig2)
        colormap(cmap_fig1)
        axis square
        hold on
        plot( RetinaAllX(ind_od_onoff==3) , RetinaAllY(ind_od_onoff==3), 'r.' ,'MarkerSize', MarkerSize)
        plot( RetinaAllX(ind_od_onoff==4), RetinaAllY(ind_od_onoff==4), 'b.' ,'MarkerSize', MarkerSize)
        title('Ipsi','fontsize',20)
        xlim([col_range_box(1) col_range_box(end)]),ylim([row_range_box(1) row_range_box(end)])
        axis on 
        
        sgtitle('All the Afferents','fontsize',20)
    end 
    %% showing just the selected afferent on colormap

    retina_contra_on_x = ret_x_aff_select(retina_od==1 & retina_onoff==1);
    retina_contra_on_y = ret_y_aff_select(retina_od==1 & retina_onoff==1);
    retina_contra_off_x = ret_x_aff_select(retina_od==1 & retina_onoff==-1);
    retina_contra_off_y = ret_y_aff_select(retina_od==1 & retina_onoff==-1);
    retina_ipsi_on_x = ret_x_aff_select(retina_od==-1 & retina_onoff==1);
    retina_ipsi_on_y = ret_y_aff_select(retina_od==-1 & retina_onoff==1);
    retina_ipsi_off_x = ret_x_aff_select(retina_od==-1 & retina_onoff==-1);
    retina_ipsi_off_y = ret_y_aff_select(retina_od==-1 & retina_onoff==-1);
    if show_fig == 1
        MarkerSize = 20;
        figure, %imshow(Retinotopy_map_plot_fig1,cmap_fig1)
        subplot(121)
        imagesc(ret_map_fig1)
        colormap(cmap_fig1)
        axis square
        hold on
        plot( retina_contra_on_x(:) , retina_contra_on_y(:) ,'r.' ,'MarkerSize', MarkerSize)
        plot( retina_contra_off_x(:), retina_contra_off_y(:),'b.' ,'MarkerSize', MarkerSize)
        title('Contra','fontsize',20)
        %xlim([0 axis_lim]),ylim([0 axis_lim])
        axis on
        
        subplot(122)
        imagesc(ret_map_fig1)
        colormap(cmap_fig1)
        axis square
        hold on
        plot( retina_ipsi_on_x(:) , retina_ipsi_on_y(:) ,'r.' ,'MarkerSize', MarkerSize)
        plot( retina_ipsi_off_x(:), retina_ipsi_off_y(:),'b.' ,'MarkerSize', MarkerSize)
        title('Ipsi','fontsize',20)
        %xlim([0 axis_lim]),ylim([0 axis_lim])
        axis on
        
        sgtitle('Selected Afferents','fontsize',20)
    end
end 


function show_onoff_od_retinotopy(OD_ONOFFSort_input,ONOFF_ONOFFSort_input,row_range,col_range,...
    length_crtx_present,Retinotopy3mIndPlot,RetinotopyRFspace_plot,ret_map_fig1,cmap_fig1)   
%% afferents after sorting OD (Orange : Ipsi, Black: Contra)

 CrtxLength = length_crtx_present; 
% row_range = 1:8; 
% col_range = 1:8; 
CircleDiameter = 20; 
OD_ODSort_small = OD_ONOFFSort_input(row_range,col_range); 

figure
subplot(131),cla
%set(gca,'units','normalized','outerposition',[0 0 1 1])
imshow(ones(CrtxLength))
title('OD Sorting','fontsize',20)
ax1 = gca; ax1.TickDir ='out';
axis square, hold on 
CrtxOD3mSorted = OD_ODSort_small; 

for g1 = 1:CrtxLength
    for g2 = 1:CrtxLength
        odTemp = CrtxOD3mSorted(g1,g2); 
        if odTemp > 0.5 % Contra
                plot( g2 , g1 ,'o' ,'MarkerSize', CircleDiameter,'LineWidth',3,'MarkerEdgeColor',[0 0 0]/255,'MarkerFaceColor',[0 0 0]/255)
        elseif odTemp < 0.5   % Ipsi
                plot( g2 , g1 ,'o' ,'MarkerSize', CircleDiameter,'LineWidth',3,'MarkerEdgeColor',[245 134 52]/255,'MarkerFaceColor',[245 134 52]/255)
        end
    end 
end 
%%   after sorting ONOFF 
OD_ONOFFSort_small = OD_ONOFFSort_input(row_range,col_range); 
ONOFF_ONOFFSort_small = ONOFF_ONOFFSort_input(row_range,col_range); 

subplot(132),cla
%set(gca,'units','normalized','outerposition',[0 0 1 1])

imshow(ones(CrtxLength))
title('ONOFF OD Sorting','fontsize',20)
ax1 = gca; 
ax1.TickDir ='out';
axis square, hold on 
CrtxOD3mSorted = OD_ONOFFSort_small; 
CrtxONOFF3mSorted = ONOFF_ONOFFSort_small; 

for g1 = 1:CrtxLength
    for g2 = 1:CrtxLength
        polTemp = CrtxONOFF3mSorted(g1,g2); 
        odTemp = CrtxOD3mSorted(g1,g2); 
        if odTemp > 0.5 % Contra
            if polTemp > 0.5  %ON
                plot( g2 , g1 ,'o' ,'MarkerSize', CircleDiameter,'LineWidth',3,'MarkerEdgeColor',[0 0 0]/255,'MarkerFaceColor',[255 0 0]/255)
            elseif polTemp < 0.5 %OFF
                plot( g2 , g1 ,'o' ,'MarkerSize', CircleDiameter,'LineWidth',3,'MarkerEdgeColor',[0 0 0]/255,'MarkerFaceColor',[0 0 255]/255)
            end 
        elseif odTemp < 0.5   % Ipsi
            if polTemp > 0.5  % ON
                plot( g2 , g1 ,'o' ,'MarkerSize', CircleDiameter,'LineWidth',3,'MarkerEdgeColor',[245 134 52]/255,'MarkerFaceColor',[255 0 0]/255)
            elseif polTemp < 0.5  % OFF
                plot( g2 , g1 ,'o' ,'MarkerSize', CircleDiameter,'LineWidth',3,'MarkerEdgeColor',[245 134 52]/255,'MarkerFaceColor',[0 0 255]/255)
            end 
        end
    end 
end 
%%  Retinotopy before sorting 
%CircleDiameter = 30;
subplot(133),cla
imshow(ones(CrtxLength))

Retinotopy_small = Retinotopy3mIndPlot(row_range,col_range); 
title('Retinotopy Before Sorting','fontsize',20)
ax1 = gca; ax1.TickDir ='out';
axis square, hold on
for g1 = 1:CrtxLength
    for g2 = 1:CrtxLength
        iColorRetinotopy = Retinotopy_small(g1,g2);
        [rr1,cc1] = find(RetinotopyRFspace_plot == iColorRetinotopy);
        i_color = ret_map_fig1(rr1,cc1);
        plot( g2 , g1 ,'o' ,'MarkerSize', CircleDiameter,'LineWidth',3,'MarkerEdgeColor',cmap_fig1(i_color,:),'MarkerFaceColor',cmap_fig1(i_color,:))
        %plot( g2 , g1 ,'o' ,'MarkerSize', CircleDiameter,'LineWidth',0.005,'MarkerEdgeColor',cmap_fig1(i_color,:),'MarkerFaceColor',cmap_fig1(i_color,:))
        %plot( g2 , g1 ,'o' ,'MarkerSize', CircleDiameter,'LineStyle','none','MarkerFaceColor',cmap(iColorRetinotopy,:))
    end
end
%axis on
end 


%% supplmentary figure1 
function show_retinotopy_with_receptive_field(row_range,col_range,length_crtx_present,Retinotopy3mIndPlot,RetinotopyRFspace_plot,ret_map_fig1,cmap_fig1,RFONCenter)
    %%  Retinotopy before sorting with receptive fields before spread 
    CrtxLength = length_crtx_present; 
    figure, %imshow(ones(CrtxLength))
    
    row_cells = 1 : RFONCenter : (RFONCenter*CrtxLength); 
    col_cells = 1 : RFONCenter : (RFONCenter*CrtxLength); 

    Retinotopy_small = Retinotopy3mIndPlot(row_range,col_range); 
    
    for g1 = 1:CrtxLength
        for g2 = 1:CrtxLength
            iColorRetinotopy = Retinotopy_small(g1,g2);
            [rr1,cc1] = find(RetinotopyRFspace_plot == iColorRetinotopy);
            i_color = ret_map_fig1(rr1,cc1);
            %plot(gca, col_cells(g2) , row_cells(g1) ,'o' ,'MarkerSize', RFONCenter,'LineWidth',3,'MarkerEdgeColor',cmap_fig1(i_color,:),'MarkerFaceColor',[1 1 1]); hold on
            viscircles([col_cells(g2) , row_cells(g1)],RFONCenter,'color', cmap_fig1(i_color,:))
            % viscircles is not supported in app designer 
            % rectangle('Position',[col_cells(g2)-RFONCenter, row_cells(g1)-RFONCenter, 2*RFONCenter, 2*RFONCenter],'EdgeColor', cmap_fig1(i_color,:),'LineWidth', 2, 'curvature', [1 1])
        end
    end
    set(gca,'ydir','reverse')
    title('Retinotopy Before Sorting with RF','fontsize',20)
    ax1 = gca; ax1.TickDir ='out';
    axis square, hold on
    xlim([-10 50]), ylim([-10 50])
end 






