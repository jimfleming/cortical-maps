function [appdata] = sort_afferent_od2_app(appdata)
    filt_sort = appdata.SortFiltOD;
    Crtx_plt_OD= appdata.CrtxOD3m;
    Crtx_plt_ONOFF= appdata.CrtxONOFF3m;
    ret_initial= appdata.Retinotopy3mIndex;
    retinotopyPlot= appdata.Retinotopy3mIndPlot;
    N_repeat= appdata.NSortOD;
    cmap= appdata.cmap;
    rng_trial= appdata.rng_trial;
    show_fig = 0;

    rng(rng_trial)

    result_OD_sorted = Crtx_plt_OD ;
    Crtx_plt_ONOFF_OD_sorted = Crtx_plt_ONOFF;
    ret_sorted = ret_initial;

    retinotopyPlot_initial = retinotopyPlot ;
    OD_holder{N_repeat} = [];
    ONOFF_holder{N_repeat} = [];
    retinotopy_holder{N_repeat} = [];
    retinotopy_plot_holder{N_repeat} = [];
    for n = 1 : N_repeat
        %   Selection of points should be randomized
        [Affr,Affc] = find( result_OD_sorted ~= 0 );
        ind_rand2 =  randperm(length(Affr)); %1 : length(Affr);
        for q = 1 : length(Affr)
            ind_replace = [];
            i1 = Affr(ind_rand2(q));
            j1 = Affc(ind_rand2(q));
            temp = zeros(size(result_OD_sorted,1),size(result_OD_sorted,2));
            temp(i1,j1) = 1;
            temp_mask1 = filter2(filt_sort,temp);
            Conv_CC = (result_OD_sorted > 0) .* temp_mask1 .* (temp_mask1>0);  % Convolution Center Contra
            Sum_CC = sum(abs(Conv_CC(:)));
            
            Conv_CI = (result_OD_sorted < 0) .* temp_mask1 .* (temp_mask1>0);  % Convolution Center Ipsi
            Sum_CI = sum(abs(Conv_CI(:)));
            
            Conv_SC = (result_OD_sorted > 0) .* temp_mask1 .* (temp_mask1<0);    % Convolution Surround Contra
            Sum_SC = sum(abs(Conv_SC(:)));
            
            Conv_SI = (result_OD_sorted < 0) .* temp_mask1 .* (temp_mask1<0);    % Convolution Surround Ipsi
            Sum_SI = sum(abs(Conv_SI(:)));
            if result_OD_sorted(i1,j1) == 1
                Conv_Result = (Sum_CC - Sum_CI) - (Sum_SC - Sum_SI);
                mask_search = filter2(ones(3),temp) - temp;
                [srch_r,srch_c] = find( mask_search & result_OD_sorted ~= 1 ); %   figure,imagesc(mask_search & result_OD_sorted ~= 1 )
                
                Conv_temp = zeros(size(srch_r,1),1);
                Aff_to_be_replaced = zeros(size(srch_r,1),1);
                
                for k1 = 1 : size(srch_r,1)
                    temp_result = result_OD_sorted;
                    temp_result(i1,j1) = result_OD_sorted(srch_r(k1),srch_c(k1));
                    temp_result(srch_r(k1),srch_c(k1)) = 1;
                    temp_mask_Contra = circshift(temp_mask1, [srch_r(k1)-i1 srch_c(k1)-j1]);
                    
                    Conv_CC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra>0)));  % Convolution Center Contra
                    Conv_CI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra>0)));  % Convolution Center Ipsi
                    Conv_SC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra<0)));    % Convolution Surround Contra
                    Conv_SI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra<0)));    % Convolution Surround Ipsi
                    Conv_temp(k1) = (Conv_CC_temp - Conv_CI_temp) - (Conv_SC_temp - Conv_SI_temp);
                    Aff_to_be_replaced(k1) = result_OD_sorted(srch_r(k1),srch_c(k1));
                end
            elseif result_OD_sorted(i1,j1) == -1
                Conv_Result = (Sum_CI - Sum_CC) - (Sum_SI - Sum_SC);
                
                mask_search = filter2(ones(3),temp) - temp; %   figure,imagesc(mask_search & result_binary ~= -1 )
                [srch_r,srch_c] = find( mask_search & result_OD_sorted ~= -1 );
                
                Conv_temp = zeros(size(srch_r,1),1);
                Aff_to_be_replaced = zeros(size(srch_r,1),1);

                for k1 = 1 : size(srch_r,1)
                    temp_result = result_OD_sorted;
                    temp_result(i1,j1) = result_OD_sorted(srch_r(k1),srch_c(k1));
                    temp_result(srch_r(k1),srch_c(k1)) = 1;
                    temp_mask_Contra = circshift(temp_mask1, [srch_r(k1)-i1 srch_c(k1)-j1]);
                    
                    Conv_CC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra>0)));  % Convolution Center Contra
                    Conv_CI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra>0)));  % Convolution Center Ipsi
                    Conv_SC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra<0)));    % Convolution Surround Contra
                    Conv_SI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra<0)));    % Convolution Surround Ipsi
                    Conv_temp(k1) = (Conv_CI_temp - Conv_CC_temp) - (Conv_SI_temp - Conv_SC_temp);
                    Aff_to_be_replaced(k1) = result_OD_sorted(srch_r(k1),srch_c(k1));
                end
            end
            if sum( Conv_temp > Conv_Result )
                ind_replace_vec = find(max(Conv_temp) == Conv_temp);
                ind_replace = ind_replace_vec(1);
                result_OD_sorted(srch_r(ind_replace),srch_c(ind_replace)) = result_OD_sorted(i1,j1);
                result_OD_sorted(i1,j1) = Aff_to_be_replaced(ind_replace);
                
                Crtx_plt_ONOFF_OD_temp = Crtx_plt_ONOFF_OD_sorted(i1,j1);
                Crtx_plt_ONOFF_OD_sorted(i1,j1) = Crtx_plt_ONOFF_OD_sorted(srch_r(ind_replace),srch_c(ind_replace));
                Crtx_plt_ONOFF_OD_sorted(srch_r(ind_replace),srch_c(ind_replace)) = Crtx_plt_ONOFF_OD_temp;
                
                ret_temp = ret_sorted(srch_r(ind_replace),srch_c(ind_replace));
                ret_sorted(srch_r(ind_replace),srch_c(ind_replace)) = ret_sorted(i1,j1);
                ret_sorted(i1,j1) = ret_temp;
                
                
                ret_temp2 = retinotopyPlot(srch_r(ind_replace),srch_c(ind_replace));
                retinotopyPlot(srch_r(ind_replace),srch_c(ind_replace)) = retinotopyPlot(i1,j1);
                retinotopyPlot(i1,j1) = ret_temp2;
                
            end
        end
        
        OD_holder{n} =  result_OD_sorted;
        ONOFF_holder{n} =  Crtx_plt_ONOFF_OD_sorted;
        retinotopy_holder{n} =  ret_sorted;
        retinotopy_plot_holder{n} =  retinotopyPlot;
    end

    map_step_holder.OD_holder = OD_holder;
    map_step_holder.ONOFF_holder = ONOFF_holder;
    map_step_holder.retinotopy_holder = retinotopy_plot_holder;

    appdata.OD_ODSort = result_OD_sorted;
    appdata.ONOFF_ODSort = Crtx_plt_ONOFF_OD_sorted;
    appdata.RetODsorted = ret_sorted;
    appdata.RetinotopyODSortedPlot = retinotopyPlot;
    appdata.map_step_holder = map_step_holder;
end
