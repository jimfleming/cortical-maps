function [appdata] = sort_afferent_onoff2_app(appdata)
    filt_sort = appdata.SortFiltONOFF;
    onoff_input = appdata.ONOFF_ODSort;
    od_input = appdata.OD_ODSort;
    retinotopy_input = appdata.RetODsorted;
    retinotopyPlot = appdata.RetinotopyODSortedPlot;
    N_repeat = appdata.NSortONOFF;
    cmap = appdata.cmap;
    rng_trial = appdata.rng_trial;

    % The same functionality as fun_aff_sort but the ON OFF can only move
    % within the same OD band
    rng(rng_trial)

    onoff_sorted = onoff_input ;
    od_sorted = od_input;
    ret_sorted = retinotopy_input;

    retinotopyPlot_initial = retinotopyPlot ;
    OD_holder{N_repeat} = [];
    ONOFF_holder{N_repeat} = [];
    retinotopy_plot_holder{N_repeat} = [];
    retinotopy_holder{N_repeat} = [];
    for n = 1 : N_repeat
        %   Selection of points should be randomized
        [Affr,Affc] = find( onoff_sorted ~= 0 );
        ind_rand2 = randperm(length(Affr));
        for q = 1 : length(Affr)
            ind_replace = [];
            i1 = Affr(ind_rand2(q));
            j1 = Affc(ind_rand2(q));
            
            temp = zeros(size(onoff_sorted,1),size(onoff_sorted,2));
            temp(i1,j1) = 1;
            temp_mask1 = filter2(filt_sort,temp);
            Conv_CC = (onoff_sorted > 0) .* temp_mask1 .* (temp_mask1>0);  % Convolution Center Contra
            Sum_CC = sum(abs(Conv_CC(:)));
            
            Conv_CI = (onoff_sorted < 0) .* temp_mask1 .* (temp_mask1>0);  % Convolution Center Ipsi
            Sum_CI = sum(abs(Conv_CI(:)));
            
            Conv_SC = (onoff_sorted > 0) .* temp_mask1 .* (temp_mask1<0);    % Convolution Surround Contra
            Sum_SC = sum(abs(Conv_SC(:)));
            
            Conv_SI = (onoff_sorted < 0) .* temp_mask1 .* (temp_mask1<0);    % Convolution Surround Ipsi
            Sum_SI = sum(abs(Conv_SI(:)));
            
            if onoff_sorted(i1,j1) == 1
                Conv_Result = (Sum_CC - Sum_CI) - (Sum_SC - Sum_SI);
                mask_search = filter2(ones(3),temp) - temp;
                if od_sorted(i1,j1) == 1
                    [srch_r,srch_c] = find( mask_search & onoff_sorted ~= 1 & od_sorted == 1); %   figure,imagesc(mask_search & result_OD_sorted ~= 1 )
                elseif od_sorted(i1,j1) == -1
                    [srch_r,srch_c] = find( mask_search & onoff_sorted ~= 1 & od_sorted == -1); %   figure,imagesc(mask_search & result_OD_sorted ~= 1 )
                end
                Conv_temp = zeros(size(srch_r,1),1);
                Aff_to_be_replaced = zeros(size(srch_r,1),1);
                for k1 = 1 : size(srch_r,1)
                    temp_result = onoff_sorted;
                    temp_result(i1,j1) = onoff_sorted(srch_r(k1),srch_c(k1));
                    temp_result(srch_r(k1),srch_c(k1)) = 1;
                    temp_mask_Contra = circshift(temp_mask1, [srch_r(k1)-i1 srch_c(k1)-j1]);
                    
                    Conv_CC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra>0)));  % Convolution Center Contra
                    Conv_CI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra>0)));  % Convolution Center Ipsi
                    Conv_SC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra<0)));    % Convolution Surround Contra
                    Conv_SI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra<0)));    % Convolution Surround Ipsi
                    Conv_temp(k1) = (Conv_CC_temp - Conv_CI_temp) - (Conv_SC_temp - Conv_SI_temp);
                    Aff_to_be_replaced(k1) = onoff_sorted(srch_r(k1),srch_c(k1));
                end
                
            elseif onoff_sorted(i1,j1) == -1
                Conv_Result = (Sum_CI - Sum_CC) - (Sum_SI - Sum_SC);
                mask_search = filter2(ones(3),temp) - temp; %   figure,imagesc(mask_search & result_binary ~= -1 )
                if od_sorted(i1,j1) == 1
                    [srch_r,srch_c] = find( mask_search & onoff_sorted ~= -1 & od_sorted == 1);
                elseif od_sorted(i1,j1) == -1
                    [srch_r,srch_c] = find( mask_search & onoff_sorted ~= -1 & od_sorted == -1);
                end
                Conv_temp = zeros(size(srch_r,1),1);
                Aff_to_be_replaced = zeros(size(srch_r,1),1);
                for k1 = 1 : size(srch_r,1)
                    temp_result = onoff_sorted;
                    temp_result(i1,j1) = onoff_sorted(srch_r(k1),srch_c(k1));
                    temp_result(srch_r(k1),srch_c(k1)) = 1;
                    temp_mask_Contra = circshift(temp_mask1, [srch_r(k1)-i1 srch_c(k1)-j1]);
                    
                    Conv_CC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra>0)));  % Convolution Center Contra
                    Conv_CI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra>0)));  % Convolution Center Ipsi
                    Conv_SC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra<0)));    % Convolution Surround Contra
                    Conv_SI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra<0)));    % Convolution Surround Ipsi
                    Conv_temp(k1) = (Conv_CI_temp - Conv_CC_temp) - (Conv_SI_temp - Conv_SC_temp);
                    Aff_to_be_replaced(k1) = onoff_sorted(srch_r(k1),srch_c(k1));
                end
            end
            
            if sum( Conv_temp > Conv_Result )
                ind_replace_vec = find(max(Conv_temp) == Conv_temp);
                ind_replace = ind_replace_vec(1);
                onoff_sorted(srch_r(ind_replace),srch_c(ind_replace)) = onoff_sorted(i1,j1);
                onoff_sorted(i1,j1) = Aff_to_be_replaced(ind_replace);
                
                Crtx_plt_ONOFF_OD_temp = od_sorted(i1,j1);
                od_sorted(i1,j1) = od_sorted(srch_r(ind_replace),srch_c(ind_replace));
                od_sorted(srch_r(ind_replace),srch_c(ind_replace)) = Crtx_plt_ONOFF_OD_temp;
                
                ret_temp = ret_sorted(srch_r(ind_replace),srch_c(ind_replace));
                ret_sorted(srch_r(ind_replace),srch_c(ind_replace)) = ret_sorted(i1,j1);
                ret_sorted(i1,j1) = ret_temp;
                
                ret_temp2 = retinotopyPlot(srch_r(ind_replace),srch_c(ind_replace));
                retinotopyPlot(srch_r(ind_replace),srch_c(ind_replace)) = retinotopyPlot(i1,j1);
                retinotopyPlot(i1,j1) = ret_temp2;
                
            end
        end
        
        OD_holder{n} =  od_sorted;
        ONOFF_holder{n} =  onoff_sorted;
        retinotopy_holder{n} =  ret_sorted;
        retinotopy_plot_holder{n} =  retinotopyPlot;
    end

    map_step_holder.OD_holder = OD_holder;
    map_step_holder.ONOFF_holder = ONOFF_holder;
    map_step_holder.retinotopy_holder = retinotopy_plot_holder;

    appdata.ONOFFCrtxPlt = onoff_sorted;
    appdata.ODCrtxPlt = od_sorted;
    appdata.RetONOFFsorted = ret_sorted;
    appdata.RetinotopyONOFFSortedPlot = retinotopyPlot;
    appdata.ONOFFmap_step_holder = map_step_holder;
end
