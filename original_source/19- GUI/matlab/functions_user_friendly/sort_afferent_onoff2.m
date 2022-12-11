function     [onoff_sorted,od_sorted,ret_sorted,retinotopyPlot,map_step_holder] = sort_afferent_onoff2(filt_sort,onoff_input,od_input,retinotopy_input,retinotopyPlot,N_repeat,cmap,rng_trial,show_fig)
    % The same functionality as fun_aff_sort but the ON OFF can only move 
    % within the same OD band 
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            %     figure, imagesc(temp_mask1)
            [Aff_to_be_replaced, Conv_temp, Conv_Result, srch_r, srch_c] = get_aff_to_be_replaced(temp, temp_mask1, od_sorted, onoff_sorted);
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


    %      Debug : OD ONOFF After Sorting ONOFF 
    if show_fig == 1 
        figure
        set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
        subplot(131),imagesc(od_sorted), colormap(gca, 'gray'), axis square,title('Ocular Dominance')
        subplot(132),imagesc(onoff_sorted), colormap(gca, 'jet'), axis square,title('ON OFF')
        subplot(133),imshow(retinotopyPlot,cmap),title('Retinotopy')
        sgtitle('ON/OFF sorting ','fontsize',30)
        
        %% sorting process 
%         figure 
%         set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
%         subplot(3,11,1),imagesc(od_input), colormap(gca,'gray'), axis square, ylabel('OD'), set(gca,'box','off','XTick',[], 'YTick', [])
%         subplot(3,11,12),imagesc(onoff_input), colormap(gca,'jet'), axis square, ylabel('ON OFF'),  set(gca,'box','off','XTick',[], 'YTick', [])
%         subplot(3,11,23),imshow(retinotopyPlot_initial,cmap), ylabel('Retinotopy'),  set(gca,'box','off','XTick',[], 'YTick', [])
%         for nn = 1 : N_repeat
%             subplot(3,11,1+nn),imagesc(OD_holder{nn}), colormap(gca, 'gray'), axis square , axis off 
%             title(sprintf('step %.0f',nn))
%             subplot(3,11,12+nn),imagesc(ONOFF_holder{nn}), colormap(gca, 'jet'), axis square, axis off 
%             subplot(3,11,23+nn),imshow(retinotopy_plot_holder{nn},cmap), axis off 
%         end
%         sgtitle('ON/OFF sorting ','fontsize',30)
        
        %% suppplementary figure 
%         cs = 10; % crop_size
%         rr = 5 : 5+cs ; 
%         cc = 5 : 5+cs ; 
%         ind_aff = 123 ;% ret_initial(1); 
%         [rrr,ccc] = find(ind_aff == retinotopy_input); 
%         figure 
%         set(gcf, 'units', 'normalized','OuterPosition',[0 0 .9 .9]);
%         subplot(1,11,1),imagesc(onoff_input(rr,cc)), colormap(gca,'jet'), axis square, ylabel('OD'), set(gca,'box','off','XTick',[], 'YTick', [])
%         hold on, plot(ccc,rrr,'ro')
%         for nn = 1 : N_repeat
%             subplot(1,11,1+nn),imagesc(ONOFF_holder{nn}(rr,cc)), colormap(gca,'jet'), axis square , axis off
%             title(sprintf('step %.0f',nn))
%             [rrr,ccc] = find(ind_aff == retinotopy_holder{nn});
%             hold on, plot(ccc,rrr,'ro')
%         end
%         sgtitle('ONOFF sorting ','fontsize',30)
%         annotation('textbox',[0.1 0.88 0.1 0.1],'string','red : on, blue : off','edgecolor','none','fontsize',14)

        
    end

end 