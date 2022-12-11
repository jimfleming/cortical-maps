function [weight_map] = synaptic_weight_island2(ori_map,ONOFFODLabelSorted,ODCrtxPlt,ONOFFCrtxPlt,show_fig)
    %   Generating synaptic weight besed on distance to pinwheels and to the
    %   border of the islands 

    %%
    weight_border = 1; %  maximum weight 
    weight_pw = 0.5;  % %  minimum weight 
    weight_map = zeros(size(ODCrtxPlt));

    %% finding the location of pinwheels 
    [pw_locations] = (get_pinwheel_location_and_charge_coverage(ori_map,0,0)); 

    %%
    od_edge = edge(ODCrtxPlt>0.5,'canny'); 
    onoff_edge = edge(ONOFFCrtxPlt>0.5,'canny'); 

    island_edge = od_edge | onoff_edge; 
    [row_border,col_border] = find(island_edge);

    for i_row = 1 : size(ODCrtxPlt,1)
        for i_col = 1 : size(ODCrtxPlt,2)
            [dist_pw,i_pw] = min(sqrt((i_row - pw_locations(:,2)).^2 +  (i_col - pw_locations(:,1)).^2 ));
            [dist_border,~] = min(sqrt((i_row - row_border(:)).^2 +  (i_col - col_border(:)).^2 ));
            weight_map(i_row,i_col) = (weight_pw * dist_border + weight_border * dist_pw) / (dist_border + dist_pw); 
        end 
    end 

    %% Assigning weight to unassigned locations 
    weight_map_complete = weight_map; 
    [r4,c4] = find( island_edge ); %& SelectedArea);
    for pp = 1 : length(r4)
        r5 = r4(pp);
        c5 = c4(pp);
        Temp3= zeros(size(ONOFFODLabelSorted));
        Temp3(r5,c5) = 1;
        Temp3 = imdilate(Temp3,ones(3)) - Temp3;
        value_available = weight_map_complete(Temp3 == 1 & ~island_edge);
        if isempty(value_available)
            continue
        else
            ave_neighbior = mean(value_available);
            weight_map_complete(r5,c5) = ave_neighbior;
        end
    end 

    weight_map = imgaussfilt(weight_map_complete, 1); 
    weight_map = weight_map - min(weight_map(:)); 
    weight_map = weight_map / max(weight_map(:)); 

    %%
    if show_fig ==  1 
        figure; imagesc(weight_map); colormap(gca, 'jet'); axis square; hold on;colorbar
        od_contour_levels = 1; 
        hold on, contour(ODCrtxPlt, od_contour_levels, 'k', 'LineWidth', 5);
        title('synaptic weight based on islands','fontsize',20)
    end 
end 
