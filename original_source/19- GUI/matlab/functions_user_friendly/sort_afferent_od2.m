function [result_OD_sorted,Crtx_plt_ONOFF_OD_sorted,ret_sorted,retinotopyPlot,map_step_holder] = sort_afferent_od2(filt_sort,Crtx_plt_OD,Crtx_plt_ONOFF,ret_initial,retinotopyPlot,N_repeat,cmap,rng_trial,show_fig)

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
            %     figure, imagesc(temp_mask1)
            [Aff_to_be_replaced, Conv_temp, Conv_Result, srch_r, srch_c] = get_aff_to_be_replaced(temp, temp_mask1, result_OD_sorted);
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
    
    if show_fig == 1 
       figure 
        set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
        subplot(131),imagesc(result_OD_sorted), colormap(gca, 'gray'), axis square ,title('Ocular Dominance')
        subplot(132),imagesc(Crtx_plt_ONOFF_OD_sorted), colormap(gca, 'jet'), axis square,title('ON OFF')
        subplot(133),imshow(retinotopyPlot,cmap),title('Retinotopy')
        sgtitle('OD sorting ','fontsize',30)
        
        %% sorting process 
%         figure 
%         set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
%         subplot(3,11,1),imagesc(Crtx_plt_OD), colormap(gca,'gray'), axis square, ylabel('OD'), set(gca,'box','off','XTick',[], 'YTick', [])
%         subplot(3,11,12),imagesc(Crtx_plt_ONOFF), colormap(gca,'jet'), axis square, ylabel('ON OFF'),  set(gca,'box','off','XTick',[], 'YTick', [])
%         subplot(3,11,23),imshow(retinotopyPlot_initial,cmap), ylabel('Retinotopy'),  set(gca,'box','off','XTick',[], 'YTick', [])
%         for nn = 1 : N_repeat
%             subplot(3,11,1+nn),imagesc(OD_holder{nn}), colormap(gca, 'gray'), axis square , axis off 
%             title(sprintf('step %.0f',nn))
%             subplot(3,11,12+nn),imagesc(ONOFF_holder{nn}), colormap(gca, 'jet'), axis square, axis off 
%             subplot(3,11,23+nn),imshow(retinotopy_plot_holder{nn},cmap), axis off 
%         end
%         sgtitle('OD sorting ','fontsize',30)
        
        %% suppplementary figure 
%         cs = 10; % crop_size
%         ind_aff = 2 ;% ret_initial(1); 
%         [rrr,ccc] = find(ind_aff == ret_initial); 
%         figure 
%         set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
%         subplot(1,11,1),imagesc(Crtx_plt_OD(1:cs,1:cs)), colormap(gca,[0 0 0;255 140 0]/255), axis square, ylabel('OD'), set(gca,'box','off','XTick',[], 'YTick', [])
%         hold on, plot(ccc,rrr,'ro')
%         for nn = 1 : N_repeat
%             subplot(1,11,1+nn),imagesc(OD_holder{nn}(1:cs,1:cs)), colormap(gca, [0 0 0;255 140 0]/255), axis square , axis off
%             title(sprintf('step %.0f',nn))
%             [rrr,ccc] = find(ind_aff == retinotopy_holder{nn});
%             hold on, plot(ccc,rrr,'ro')
%         end
%         sgtitle('OD sorting ','fontsize',30)
%         annotation('textbox',[0.1 0.88 0.1 0.1],'string','orange : contra, black : ipsi','edgecolor','none','fontsize',14)
    end 
end 
 
