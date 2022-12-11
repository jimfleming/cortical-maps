% supplementary figure 1, showing the sorting process

filt_sort = functionDoGfilter(2,2,5,5,1);
show_fig = 1;

%% od 
[XX, YY] = meshgrid(1:60);
od_initial = (XX < 31) * -1 + (XX >= 31) * 1;
od_initial(30, 31) = -1 ;
od_initial(30, 30) = 1 ;


%sort_afferents(od_initial, filt_sort, 'gray', show_fig)
colormap = [0 0 0; 254 134 52] / 255 ; % black and orange 
sort_afferents(od_initial, filt_sort, colormap, show_fig)
%% onoff 
onoff_initial = (YY < 31) * -1 + (YY >= 31) * 1;
onoff_initial(31, 30) = -1 ;
onoff_initial(30, 30) = 1 ;

sort_afferents(onoff_initial, filt_sort, 'jet', show_fig)


%%

function sort_afferents(input, filt_sort, cmap, show_fig)
result_OD_sorted = input ;

Affr = 30 ;
Affc = 30 ;

[srch_c,srch_r] = meshgrid(-1:1);

ind_replace = [];
temp = zeros(size(result_OD_sorted,1),size(result_OD_sorted,2));
temp(Affr,Affc) = 1;
temp_mask1 = filter2(filt_sort,temp);

%     Conv_CC = (result_OD_sorted > 0) .* temp_mask1 .* (temp_mask1>0);    % Convolution Center Contra
%     Sum_CC = sum(abs(Conv_CC(:)));
%
%     Conv_CI = (result_OD_sorted < 0) .* temp_mask1 .* (temp_mask1>0);    % Convolution Center Ipsi
%     Sum_CI = sum(abs(Conv_CI(:)));
%
%     Conv_SC = (result_OD_sorted > 0) .* temp_mask1 .* (temp_mask1<0);    % Convolution Surround Contra
%     Sum_SC = sum(abs(Conv_SC(:)));
%
%     Conv_SI = (result_OD_sorted < 0) .* temp_mask1 .* (temp_mask1<0);    % Convolution Surround Ipsi
%     Sum_SI = sum(abs(Conv_SI(:)));

if result_OD_sorted(Affr,Affc) == 1
    %        Conv_Result = (Sum_CC - Sum_CI) - (Sum_SC - Sum_SI);
    Conv_Result = sum(result_OD_sorted .* temp_mask1, [1, 2]) ;
    %         mask_search = filter2(ones(3),temp) - temp;
    %         [srch_r,srch_c] = find( mask_search & result_OD_sorted ~= 1 ); %   figure,imagesc(mask_search & result_OD_sorted ~= 1 )
    
    Conv_temp = zeros(size(srch_r,1),1);
    Aff_to_be_replaced = zeros(size(srch_r,1),1);
    
    for k1 = 1 : length(srch_r(:))
        temp_result = result_OD_sorted;
        temp_result(Affr,Affc) = result_OD_sorted(Affr + srch_r(k1), Affc + srch_c(k1));
        temp_result(Affr + srch_r(k1), Affc + srch_c(k1)) = 1;
        temp_mask_Contra = circshift(temp_mask1, [  srch_r(k1),   srch_c(k1)]);
        
        %             Conv_CC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra>0)));    % Convolution Center Contra
        %             Conv_CI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra>0)));    % Convolution Center Ipsi
        %             Conv_SC_temp = sum(sum((temp_result > 0) .* temp_mask_Contra .* (temp_mask_Contra<0)));    % Convolution Surround Contra
        %             Conv_SI_temp = sum(sum((temp_result < 0) .* temp_mask_Contra .* (temp_mask_Contra<0)));    % Convolution Surround Ipsi
        %             Conv_temp(k1) = (Conv_CC_temp - Conv_CI_temp) - (Conv_SC_temp - Conv_SI_temp);
        
        Conv_main = sum( temp_result .* temp_mask_Contra, [1, 2]);    % Convolution Center Contra
        Conv_temp(k1) = (Conv_main ) ;
        
        Aff_to_be_replaced(k1) = result_OD_sorted(Affr + srch_r(k1), Affc + srch_c(k1));
        result_holder{k1} = temp_result ;
    end
end

if sum( Conv_temp > Conv_Result )
    ind_replace_vec = find(max(Conv_temp) == Conv_temp);
    ind_replace = ind_replace_vec(1);
    result_OD_sorted(Affr + srch_r(ind_replace), Affc + srch_c(ind_replace)) = result_OD_sorted(Affr,Affc);
    result_OD_sorted(Affr,Affc) = Aff_to_be_replaced(ind_replace);
end


%%
if show_fig == 1
    %     figure
    %     set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
    %     subplot(121), imagesc(od_initial), colormap(gca, 'gray'), axis square ,title('Ocular Dominance')
    %     subplot(122), imagesc(result_OD_sorted), colormap(gca, 'gray'), axis square ,title('Ocular Dominance')
    %     sgtitle('OD sorting ','fontsize',30)
    
    figure,
    set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
    k_plot = reshape(1:9,[3,3])';
    for k1 = 1 : length(srch_r(:))
        subplot(3,3,k_plot(k1))
        imagesc(result_holder{k1}), colormap(gca, cmap), axis square, axis off
        rectangle('Position',[Affc + srch_c(k1) - .5, Affr + srch_r(k1) - .5, 1 , 1],'EdgeColor','r', 'linewidth', 3)
        xlim([27.5 32.5]), ylim([27.5 32.5])
        title(Conv_temp(k1))
    end
end

end 




