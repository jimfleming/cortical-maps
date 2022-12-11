classdef gui_JMA_112221a < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        Panel                         matlab.ui.container.Panel
        GridLayout18                  matlab.ui.container.GridLayout
        UIAxes_ASF                    matlab.ui.control.UIAxes
        UIAxes_aff_RF_vis_space       matlab.ui.control.UIAxes
        UIAxes_RGC_ipsi               matlab.ui.control.UIAxes
        UIAxes_RGC_contra             matlab.ui.control.UIAxes
        TabGroup2                     matlab.ui.container.TabGroup
        CorticaldevelopmentTab        matlab.ui.container.Tab
        GridLayout6                   matlab.ui.container.GridLayout
        title_retinotopy              matlab.ui.control.EditField
        title_eye_polarity            matlab.ui.control.EditField
        title_segregation_od          matlab.ui.control.EditField
        title_mature_orimap           matlab.ui.control.EditField
        title_primord_orimap          matlab.ui.control.EditField
        title_eye_pol                 matlab.ui.control.EditField
        UIAxes_primord_orimap_contra  matlab.ui.control.UIAxes
        UIAxes_ONOFF_ipsi             matlab.ui.control.UIAxes
        UIAxes_ONOFF_contra           matlab.ui.control.UIAxes
        UIAxes_OD_ONOFF               matlab.ui.control.UIAxes
        UIAxes_Aff_retinotopy         matlab.ui.control.UIAxes
        UIAxes_Segregation_OD         matlab.ui.control.UIAxes
        UIAxes_mature_orimap_ipsi     matlab.ui.control.UIAxes
        UIAxes_primord_orimap_ipsi    matlab.ui.control.UIAxes
        UIAxes_mature_orimap_contra   matlab.ui.control.UIAxes
        MaturevisualcortexTab         matlab.ui.container.Tab
        GridLayout10                  matlab.ui.container.GridLayout
        UIAxes_LPI                    matlab.ui.control.UIAxes
        UIAxes_LHI                    matlab.ui.control.UIAxes
        UIAxes_retinotopy             matlab.ui.control.UIAxes
        UIAxes_SF50                   matlab.ui.control.UIAxes
        UIAxes_CV                     matlab.ui.control.UIAxes
        UIAxes3_orimap_mature         matlab.ui.control.UIAxes
        UIAxes2_ONOFF_mature          matlab.ui.control.UIAxes
        UIAxes_OD_mature              matlab.ui.control.UIAxes
        Mapdimensionrelations1Tab     matlab.ui.container.Tab
        GridLayout15                  matlab.ui.container.GridLayout
        UIAxes_slope2                 matlab.ui.control.UIAxes
        UIAxes_slope1                 matlab.ui.control.UIAxes
        UIAxes_stat_cv_lhi            matlab.ui.control.UIAxes
        UIAxes_stat_lpi_lhi           matlab.ui.control.UIAxes
        UIAxes_stat_sf_lhi            matlab.ui.control.UIAxes
        UIAxes_stat_lpi_cv            matlab.ui.control.UIAxes
        UIAxes_stat_sf_cv             matlab.ui.control.UIAxes
        Mapdimensionrelations2Tab_2   matlab.ui.container.Tab
        map_dim_relations2            matlab.ui.container.Panel
        AfferentconvergenceTab_2      matlab.ui.container.Tab
        Panel_aff_convergence         matlab.ui.container.Panel
        AfferentaxonarborsTab         matlab.ui.container.Tab
        Panel_aff_axon_arbors         matlab.ui.container.Panel
        VisualdeprivationTab          matlab.ui.container.Tab
        GridLayout16                  matlab.ui.container.GridLayout
        GridLayout_ori_biased         matlab.ui.container.GridLayout
        BiaspercentDropDown           matlab.ui.control.DropDown
        BiaspercentDropDownLabel      matlab.ui.control.Label
        GenerateButton                matlab.ui.control.Button
        OrientationbiasedEditField    matlab.ui.control.NumericEditField
        OrientationbiasedEditFieldLabel  matlab.ui.control.Label
        UIAxes_ori_od_polar_hist      matlab.ui.control.UIAxes
        UIAxes_od_biased              matlab.ui.control.UIAxes
        UIAxes_orimap_biased          matlab.ui.control.UIAxes
        UIAxes_ori_od_hist            matlab.ui.control.UIAxes
        UIAxes_orimap_contour         matlab.ui.control.UIAxes
        GridLayout14                  matlab.ui.container.GridLayout
        Y2EditFieldLabel              matlab.ui.control.Label
        Y2_rf                         matlab.ui.control.EditField
        Y1EditFieldLabel              matlab.ui.control.Label
        Y1_rf                         matlab.ui.control.EditField
        X2EditFieldLabel              matlab.ui.control.Label
        X2_rf                         matlab.ui.control.EditField
        X1EditFieldLabel              matlab.ui.control.Label
        X1_rf                         matlab.ui.control.EditField
        GridLayout11                  matlab.ui.container.GridLayout
        GridLayout12                  matlab.ui.container.GridLayout
        GridLayout13                  matlab.ui.container.GridLayout
        Lamp                          matlab.ui.control.Lamp
        Lamp_2                        matlab.ui.control.Lamp
        EyeofpreferenceSwitchLabel    matlab.ui.control.Label
        EyeofpreferenceSwitch         matlab.ui.control.ToggleSwitch
        TextArea                      matlab.ui.control.TextArea
        DrawalineonthemapButton       matlab.ui.control.Button
        UIAxes_rf_select_map          matlab.ui.control.UIAxes
        ElectrodePanel                matlab.ui.container.Panel
        GridLayout8                   matlab.ui.container.GridLayout
        TabGroup                      matlab.ui.container.TabGroup
        MainTab                       matlab.ui.container.Tab
        GridLayout3                   matlab.ui.container.GridLayout
        TextArea_processing           matlab.ui.control.TextArea
        CortexLengthDropDown          matlab.ui.control.DropDown
        CloseButton                   matlab.ui.control.Button
        StartButton                   matlab.ui.control.Button
        RetinathalamusTab             matlab.ui.container.Tab
        GridLayout20                  matlab.ui.container.GridLayout
        JitterseedDropDown            matlab.ui.control.DropDown
        JitterseedDropDownLabel       matlab.ui.control.Label
        RGCdistyDropDown              matlab.ui.control.DropDown
        RGCdistyDropDownLabel         matlab.ui.control.Label
        RGCdistxDropDown              matlab.ui.control.DropDown
        RGCdistxDropDownLabel         matlab.ui.control.Label
        AffdenfactorDropDown          matlab.ui.control.DropDown
        AffdenfactorDropDownLabel     matlab.ui.control.Label
        JitterOFFyDropDown            matlab.ui.control.DropDown
        JitterOFFyDropDownLabel       matlab.ui.control.Label
        JitterOFFxDropDown            matlab.ui.control.DropDown
        JitterOFFxDropDownLabel       matlab.ui.control.Label
        JitterONyDropDown             matlab.ui.control.DropDown
        JitterONyDropDownLabel        matlab.ui.control.Label
        JitterONxDropDown             matlab.ui.control.DropDown
        JitterONxDropDownLabel        matlab.ui.control.Label
    end


    properties (Access = public)
        % MarkerSize = 20;
        
        od_contour_w = 1; %OD  contour linewidth
        Switchvalue  = 1; %1: contra, 2: Ipsi

        %         FirstPoint  = [1 1];
        %         SecondPoint = [30 30];
        %         CurrentPoint = 1;
        %         FinalPrecessDone = false;

        appdata;
    end
    % Common functions
    methods (Access = private)


        function updatRetinaUI(app, appdata, eye_input)
            %             imshow(ones(size(appdata.rf_space_x,1)), 'Parent',axes); hold(axes,'on');
            %             plot(axes, appdata.RetONxContra(:), appdata.RetONyContra(:),'r.','MarkerSize', app.MarkerSize);
            %             plot(axes, appdata.RetOFFxContra(:), appdata.RetOFFyContra(:),'b.','MarkerSize', app.MarkerSize);
            %             xlim(axes, [2.5 32.5]), ylim(axes, [2.5 32.5])


            RetinaAllOD = appdata.RetinaAllOD;
            RetinaAllONOFF = appdata.RetinaAllONOFF;
            RetinaAllX = appdata.RetinaAllX;
            RetinaAllY = appdata.RetinaAllY;
            rf_space_x = appdata.rf_space_x;

            if strcmp(eye_input, 'contra')
                ax = app.UIAxes_RGC_contra;
                sign_od = 1;
                color_on  = [1 0 0];
                color_off = [0 0 1];
                contra_font_color = appdata.contra_font_color; % [0 0 0];
            elseif strcmp(eye_input, 'ipsi')
                ax = app.UIAxes_RGC_ipsi;
                sign_od = -1;
                color_on  = [255 100 100]/255;
                color_off = [0 175 255]/255;
                ipsi_font_color = appdata.ipsi_font_color; % [245 134 52]/255;
            end

            RetONx  = RetinaAllX(RetinaAllOD==sign_od & RetinaAllONOFF==1);
            RetONy  = RetinaAllY(RetinaAllOD==sign_od & RetinaAllONOFF==1);
            RetOFFx = RetinaAllX(RetinaAllOD==sign_od & RetinaAllONOFF==-1);
            RetOFFy = RetinaAllY(RetinaAllOD==sign_od & RetinaAllONOFF==-1);

            ax.Visible = 'on'; cla(ax)
            MarkerSize = 15 ;% 20;
            row_adj = -2.4 ;
            col_adj = -2.3 ;
            imshow(ones(size(rf_space_x,1)), 'Parent',ax)
            hold(ax, 'on'), plot(ax, RetONx(:) , RetONy(:) , '.', 'color', color_on, 'MarkerSize', MarkerSize)
            setAxes(app, ax);
            hold(ax, 'on'), plot(ax, RetOFFx(:), RetOFFy(:), '.', 'color', color_off,'MarkerSize', MarkerSize)
            setAxes(app, ax);

            %             xlim(ax, [col_min+col_adj col_min+col_adj+ret_box_length])
            %             ylim(ax, [row_min+row_adj row_min+row_adj+ret_box_length])
            xlim(ax, [2.5 32.5]), ylim(ax, [2.5 32.5])
            axis(ax, 'on'),

            % title(ax, eye_input,'fontsize',20)
            if strcmp(eye_input, 'contra')
                ax.XColor = contra_font_color;
                ax.YColor = contra_font_color;
                title(ax, sprintf('Contralateral \n retina (contra)'))
            elseif strcmp(eye_input, 'ipsi')
                ax.XColor = ipsi_font_color;
                ax.YColor = ipsi_font_color;
                title(ax, sprintf('Ipsilateral \n retina (ipsi)'), 'Color', ipsi_font_color)
            end
            pause(0.001)
        end

        function update_thalamic_aff_UI(app, appdata)
            % update thalamic afferent sampling density
            ax = app.UIAxes_aff_RF_vis_space;

            CrtxONOFF3m = appdata.CrtxONOFF3m;
            Retinotopy3mIndex = appdata.Retinotopy3mIndex;
            rf_space_x = appdata.rf_space_x;
            imshow(ones(size(rf_space_x,1)), 'Parent', ax)
            hold(ax, 'on')

            sum_aff = zeros(size(rf_space_x,1)); 
            disk_element = strel('disk', round(appdata.RFONCenter), 0);
            
            w_cortex = size(CrtxONOFF3m, 1);
            hh = round(w_cortex/2);
            range = hh - 4 : hh + 3; % 8 cortical pixels (64 afferents in total)

            row_aff_all = zeros(64, 1);
            col_aff_all = zeros(64, 1);
            cc = 0;
            for ii = range
                for jj = range
                    indRf = Retinotopy3mIndex(ii, jj);
                    [row_aff, col_aff] = find(indRf == appdata.IndRetinaRfspace);

                    if CrtxONOFF3m(ii, jj) > 0
                        appviscircles(ax, [col_aff, row_aff], appdata.RFONCenter, 'Color', 'r', 'Linewidth', 1);
                    elseif CrtxONOFF3m(ii, jj) < 0
                        appviscircles(ax, [col_aff, row_aff], appdata.RFONCenter, 'Color', 'b', 'Linewidth', 1);
                    end

                    cc = cc + 1;
                    row_aff_all(cc) = row_aff;
                    col_aff_all(cc) = col_aff;
                    
                     % to calculate the number of overlapping receptive fields in visual space 
                    temp_sum_aff = zeros(size(rf_space_x, 1)); 
                    temp_sum_aff(row_aff, col_aff) = 1; 
                    sum_aff = imdilate(temp_sum_aff, disk_element) + sum_aff;  
                end
            end
%             min_rows = min(row_aff_all(:)); 
%             max_cols  = max(col_aff_all(:)); 
            xlim(ax, [min(col_aff_all(:))-10 max(col_aff_all(:))+10])
            ylim(ax, [min(row_aff_all(:))-10 max(row_aff_all(:))+10])

            max_number_overlap_rf = max(sum_aff(:)); 
            title(ax, sprintf('Afferent sampling \n density (%.0f RF per point)', max_number_overlap_rf))

            setAxes(app, ax);
            % axis(ax, 'square')
            ax.XColor = [0 0 0];
            ax.YColor = [0 0 0];
            pause(0.001)
        end

        function [appdata] = show_maps_segregation(app, appdata)

            ONOFF_smoothed   = appdata.ONOFF_smoothed;
            ODCrtxPlt_smooth = appdata.ODCrtxPlt_smooth;
            map_retinotopy   = appdata.RetinotopyONOFFSortedPlot;
            cmap_map_retinotopy = app.appdata.cmap;
            % ONOFF_interpolated = appdata.ONOFF_interpolated;
            % ODCrtxPlt_interpolated = appdata.ODCrtxPlt_interpolated;

            % map1 = [133 94 82;255 153 153;33 31 133;153 153 255]./255;
            onoff_od_colormap = [255 100 100; 255 0 0; 0 175 255; 0 0 255]./255;
            appdata.map1 = onoff_od_colormap;
            % appdata.onoff_od_colormap = onoff_od_colormap;

            width_cortex = size(ODCrtxPlt_smooth, 1);
            %% figure
            axes_od = app.UIAxes_Segregation_OD ;
            axes_od_onoff = app.UIAxes_OD_ONOFF ;
            axes_retinotopy = app.UIAxes_Aff_retinotopy ;

            % axes_od.Visible = 'on';
            cla(axes_od), imagesc(ODCrtxPlt_smooth>.5, 'Parent', axes_od)
            hold(axes_od, 'on'),
            contour(axes_od, ODCrtxPlt_smooth, app.appdata.od_contour_levels, 'k', 'LineWidth', 1);
            caxis(axes_od, [0 1])
            colormap(axes_od, 'gray'), axis(axes_od, 'square')
            %title(axes_od, 'Segregation OD') % , 'fontsize', 14)
            app.title_segregation_od.Visible = 'on';
            xlim(axes_od, [1 width_cortex]),
            ylim(axes_od, [1 width_cortex])
            % axis(axes_od, 'off')
            setAxes(app, axes_od)

            axes_od_onoff.Visible = 'on';
            z3 = (double((ODCrtxPlt_smooth>0.5))+1)*1 ;
            z4 = (double(ONOFF_smoothed<.5)+0)*2 ;
            imagesc(z3+z4, 'Parent', axes_od_onoff)
            colormap(axes_od_onoff, onoff_od_colormap);  axis(axes_od_onoff, 'square'),
            hold(axes_od_onoff, 'on')
            % title(axes_od_onoff, 'Afferent eye-polarity') % ,'fontsize',14)
            app.title_eye_polarity.Visible = 'on';
            hold(axes_od_onoff, 'on'),
            contour(axes_od_onoff, ODCrtxPlt_smooth, app.appdata.od_contour_levels, 'k', 'LineWidth', 1);
            xlim(axes_od_onoff, [1 width_cortex]), ylim(axes_od, [1 width_cortex])
            % axis(axes_od_onoff, 'off')
            setAxes(app, axes_od_onoff)

            axes_retinotopy.Visible = 'on';
            cla(axes_retinotopy), % imagesc(map_retinotopy, 'Parent', axes_onoff)
            imshow(map_retinotopy, cmap_map_retinotopy, 'Parent', axes_retinotopy)
            % colormap(axes_onoff, cmap_map_retinotopy),
            axis(axes_retinotopy, 'square')
            %title(axes_retinotopy, 'Afferent retinotopy') % , 'fontsize', 14)
            app.title_retinotopy.Visible = 'on';
            xlim(axes_retinotopy, [1 width_cortex]), ylim(axes_retinotopy, [1 width_cortex])
            %axis(axes_onoff, 'off')
            setAxes(app, axes_retinotopy)

            pause(0.01)
        end

        function updatSmoothUI(app, appdata, axes)

            z5 = double(appdata.ODCrtxPlt_interpolated>0.5)+1 ;
            z6 = double(appdata.ONOFF_interpolated<.5)*2 ;

            imshow(z5+z6,appdata.map1, 'Parent',axes); hold(axes,'on');
            contour(axes,appdata.ODCrtxPlt_interpolated, app.appdata.od_contour_levels, 'k', 'LineWidth', 5);
            setAxes(app,axes);
        end

        function update_primord_orimap_UI(app, appdata)

            app.title_primord_orimap.Visible = 'on';

            %             app.UIAxes_primord_orimap_ipsi.Visible = 'on';
            %             axes(1) = app.UIAxes_primord_orimap_ipsi;
            %             imshow(appdata.primord_orimap_smooth, hsv(180), 'Parent',axes(1)); hold(axes(1),'on');
            %             contour(axes(1), appdata.ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
            %             setAxes(app,axes(1));

            primord_orimap_contra = appdata.data_primord_contra.OriPreferred;
            primord_orimap_ipsi = appdata.data_primord_ipsi.OriPreferred;

            axes1 = app.UIAxes_primord_orimap_contra;
            cla(axes1)
            imshow(primord_orimap_contra, hsv(180), 'Parent', axes1); hold(axes1, 'on');
            contour( axes1, appdata.ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
            % title(axes1, 'Primordial orientation map (contra)') % , 'fontsize', 14)
            xlabel(axes1, 'contra') % , 'fontsize', 14)
            setAxes(app, axes1);
            axes1.XColor = appdata.contra_font_color;
            axes1.YColor = appdata.contra_font_color;

            axes2 = app.UIAxes_primord_orimap_ipsi;
            cla(axes2)
            imshow(primord_orimap_ipsi, hsv(180), 'Parent', axes2); hold(axes2, 'on');
            contour(axes2, appdata.ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
            % title(axes2, 'Primordial orientation map (ipsi)')
            xlabel(axes2, 'ipsi')
            setAxes(app, axes2);
            axes2.XColor = appdata.ipsi_font_color;
            axes2.YColor = appdata.ipsi_font_color;

        end

        function update_mature_orimaps_UI(app, appdata)
            app.title_mature_orimap.Visible = 'on';

            mature_orimap_contra = appdata.data_contra_mature.ori_map_interpolated;
            mature_orimap_ipsi = appdata.data_ipsi_mature.ori_map_interpolated;

            ax1 = app.UIAxes_mature_orimap_contra;
            imshow(mature_orimap_contra, hsv(180), 'Parent', ax1); hold(ax1,'on');
            % title(ax1, 'Adult orientation map (contra)')
            xlabel(ax1, 'contra')
            contour(ax1, appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', 1);% app.od_contour_w);
            setAxes(app, ax1);

            ax2 = app.UIAxes_mature_orimap_ipsi;
            imshow(mature_orimap_ipsi, hsv(180), 'Parent', ax2); hold(ax2,'on');
            % title(ax2, 'Adult orientation map (ipsi)')
            xlabel(ax2, 'ipsi')
            contour(ax2, appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', 1); %app.od_contour_w);
            setAxes(app, ax2);
            ax2.XColor = appdata.ipsi_font_color;
            ax2.YColor = appdata.ipsi_font_color;
        end

        function update_mature_onoff_maps_UI(app, appdata, eye_input)

            app.title_eye_pol.Visible = 'on';
            if strcmp(eye_input, 'contra')
                % ODI = appdata.data_primord_contra.ODI;
                max_on_response = appdata.data_contra_mature.max_on_response_norm;
                max_off_response = appdata.data_contra_mature.max_off_response_norm;
                axes_onoff = app.UIAxes_ONOFF_contra;
                border_color = appdata.contra_font_color;
            elseif strcmp(eye_input, 'ipsi')
                max_on_response = appdata.data_ipsi_mature.max_on_response_norm;
                max_off_response = appdata.data_ipsi_mature.max_off_response_norm;
                axes_onoff = app.UIAxes_ONOFF_ipsi;
                border_color = appdata.ipsi_font_color;
            end

            % onoff map
            onoff_resp = max_on_response + (-1) * max_off_response;
            onoff_norm = onoff_resp / max(abs(onoff_resp(:)));
            % onoff_norm = ODI .* onoff_norm;
            onoff_norm_interp = imresize(onoff_norm, appdata.n_interpol);

            axes_onoff.Visible = 'on';
            imagesc(axes_onoff, onoff_norm_interp), axis(axes_onoff, 'square'), colormap(axes_onoff,'jet'),
            xlabel(axes_onoff, eye_input)
            hold(axes_onoff,'on')
            contour(axes_onoff, appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', 1); % app.od_contour_w);
            axis(axes_onoff, 'off')
            setAxes(app, axes_onoff);
            axes_onoff.XColor = border_color;
            axes_onoff.YColor = border_color;
        end

        function update_mature_maps_UI(app, appdata)
            n_interpol = appdata.n_interpol;

            if app.Switchvalue == 1
                ODI = appdata.data_primord_contra.ODI;

                max_on_response = appdata.data_contra_mature.max_on_response_norm;
                max_off_response = appdata.data_contra_mature.max_off_response_norm;

                Ori2 = appdata.data_contra_mature.ori_map_interpolated;
                LHI2 = appdata.data_contra_mature.LHI_interpolated;
                CV2 = appdata.data_contra_mature.cv_map_interpolated;
                LPI2 = appdata.data_contra_mature.LPI_intepolated;
                SF50_2 = appdata.data_contra_mature.sf50_map_interpolated;

            else
                ODI = 1 - appdata.data_primord_ipsi.ODI;

                max_on_response = appdata.data_ipsi_mature.max_on_response_norm;
                max_off_response = appdata.data_ipsi_mature.max_off_response_norm;

                Ori2 = appdata.data_ipsi_mature.ori_map_interpolated;
                LHI2 = appdata.data_ipsi_mature.LHI_interpolated;
                CV2 = appdata.data_ipsi_mature.cv_map_interpolated;
                LPI2 = appdata.data_ipsi_mature.LPI_intepolated;
                SF50_2 = appdata.data_ipsi_mature.sf50_map_interpolated;

            end

            map_retinotopy = appdata.RetinotopyONOFFSortedPlot;
            cmap_map_retinotopy = app.appdata.cmap;

            % retinotopy
            axes_retinotopy = app.UIAxes_retinotopy;
            axes_retinotopy.Visible = 'on';
            imshow(map_retinotopy, cmap_map_retinotopy, 'Parent', axes_retinotopy)
            title(axes_retinotopy, 'Cortical retinotopy') % , 'fontsize', 14)
            axis(axes_retinotopy, 'square'), % colormap(axes_retinotopy, 'gray'),
            axis(axes_retinotopy, 'off')
            setAxes(app, axes_retinotopy);

            % od
            ODI_interp = imresize(ODI, n_interpol);
            axes_od = app.UIAxes_OD_mature;
            axes_od.Visible = 'on';
            imagesc(axes_od, ODI_interp), axis(axes_od, 'square'), colormap(axes_od, 'gray'),
            title(axes_od, 'Ocular dominance')
            hold(axes_od,'on')
            contour(axes_od, appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', app.od_contour_w);
            axis(axes_od, 'off')
            setAxes(app, axes_od);

            % onoff map
            onoff_resp = max_on_response + (-1) * max_off_response;
            onoff_norm = onoff_resp / max(abs(onoff_resp(:)));
            %onoff_norm = ODI.* onoff_norm;
            onoff_norm_interp = imresize(onoff_norm, n_interpol);

            axes_onoff = app.UIAxes2_ONOFF_mature;
            axes_onoff.Visible = 'on';
            imagesc(axes_onoff, onoff_norm_interp), axis(axes_onoff, 'square'), colormap(axes_onoff,'jet'),
            title(axes_onoff, 'ON-OFF polarity')
            hold(axes_onoff,'on')
            contour(axes_onoff, appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', app.od_contour_w);
            axis(axes_onoff, 'off')
            setAxes(app, axes_onoff);

            %%
            NLHI2  = 255*(LHI2-min(LHI2(:)))/(max(LHI2(:))-min(LHI2(:)));
            NCV2   = 255*(CV2-min(CV2(:)))/(max(CV2(:))-min(CV2(:)));
            NLPI2  = 255*(LPI2-min(LPI2(:)))/(max(LPI2(:))-min(LPI2(:)));
            NSF50_2 = 255*(SF50_2-min(SF50_2(:)))/(max(SF50_2(:))-min(SF50_2(:)));

            axes(1) = app.UIAxes3_orimap_mature;
            imshow(Ori2, hsv(180), 'Parent',axes(1)); hold(axes(1),'on'); %freezeColors;colorbar(axes(1));
            title(axes(1), 'Orientation preference')
            contour(axes(1), appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', app.od_contour_w);
            setAxes(app,axes(1));

            axes(2) = app.UIAxes_LHI;
            imshow(NLHI2, jet(256), 'Parent',axes(2)); hold(axes(2),'on');%freezeColors;colorbar(axes(2));
            title(axes(2), 'Orientation clustering')
            contour(axes(2), appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', app.od_contour_w);
            setAxes(app,axes(2));

            axes(3) = app.UIAxes_CV;
            imshow(NCV2, jet(256), 'Parent',axes(3)); hold(axes(3),'on');%freezeColors;colorbar(axes(3));
            title(axes(3), 'Orientation selectivity')
            contour(axes(3), appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', app.od_contour_w);
            setAxes(app,axes(3));

            axes(4) = app.UIAxes_LPI;
            imshow(NLPI2, jet(256), 'Parent',axes(4)); hold(axes(4),'on');%freezeColors;colorbar(axes(4));
            title(axes(4), 'Low-pass index')
            contour(axes(4), appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', app.od_contour_w);
            setAxes(app,axes(4));

            axes(5) = app.UIAxes_SF50;
            imshow(NSF50_2, jet(256), 'Parent',axes(5)); hold(axes(5),'on');%freezeColors;colorbar(axes(5));
            title(axes(5), 'Spatial resolution')
            contour(axes(5), appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', app.od_contour_w);
            setAxes(app,axes(5));
        end

        function update_maps_intersection_UI(app, appdata)

            n_interpol = appdata.n_interpol;

            max_on_response = appdata.data_contra_mature.max_on_response_norm;
            max_off_response = appdata.data_contra_mature.max_off_response_norm;
            onoff_norm = max_on_response + (-1) * max_off_response;
            onoff_norm = onoff_norm / max(abs(onoff_norm(:)));
            onoff_interpolated = imresize(onoff_norm, n_interpol);

            [ori_map_smooth, ori_map_interpolated]  = smooth_ori_euler(appdata.data_contra_mature.ori_map_smooth, appdata.ODCrtxPlt_smooth, n_interpol, 0, 0);
            % ori_map_interpolated = appdata.data_contra_mature.ori_map_interpolated;
            SF_map_interpolated = appdata.data_contra_mature.sf50_map_interpolated;

            %% Retinotopy (Finding the center of RFs for RetinotopyX and RetinotopyY)
            row_range = 1:size(max_on_response, 1);
            col_range = 1:size(max_on_response, 2);
            [rf_center_row, rf_center_col, sigma_rfs] = find_center_sigma_rfs_S4(appdata.rf_contra_mature, row_range, col_range, 10, 25, 5, 1.75, 0);
            RetinotopyX_interpolated = imresize(rf_center_col, n_interpol);
            RetinotopyY_interpolated = imresize(rf_center_row, n_interpol);

            % map smoothing
            onoff_interpolated_smooth = imgaussfilt(onoff_interpolated,7);
            sf_interpolated_smooth = imgaussfilt(SF_map_interpolated,5);
            RetinotopyX_interpolated_smooth = imgaussfilt(RetinotopyX_interpolated,5);
            RetinotopyY_interpolated_smooth = imgaussfilt(RetinotopyY_interpolated,5);

            od_contour_interpolated = imresize(appdata.ODCrtxPlt_smooth, n_interpol);

            % t = tiledlayout(app.Mapdimensionrelations2Tab, 2, 10, 'TileSpacing', 'none', 'Padding', 'none');
            app.map_dim_relations2.Visible = 'on'; 
            cc = 11; 
            t = tiledlayout(app.map_dim_relations2, 2, cc, 'TileSpacing', 'compact', 'Padding', 'compact');
            
            % ori map
            compute_intersect_gui( ori_map_interpolated, od_contour_interpolated, 10, 1, 0, 3, hsv, 'Ori vs OD', [nexttile(t, 4) nexttile(t, 4+cc)]);
            compute_intersect_gui( ori_map_interpolated, sf_interpolated_smooth, 10, 1, 0, 3, hsv, 'Ori vs SF',  [nexttile(t, 1) nexttile(t, 1+cc)]);
            compute_intersect_gui( ori_map_interpolated, RetinotopyX_interpolated_smooth, 10, 20, 0, 3, hsv, 'Ori vs Azimuth',  [nexttile(t, 10) nexttile(t, 10+cc)]);
            compute_intersect_gui( ori_map_interpolated, RetinotopyY_interpolated_smooth, 10, 20, 0, 3, hsv, 'Ori vs Elevation',  [nexttile(t, 8) nexttile(t, 8+cc)]);

            % ONOFF (needs the caxis([-1 1])) to be added
            compute_intersect_gui(onoff_interpolated_smooth, od_contour_interpolated, 1, 1, 0, 3, jet, 'ONOFF vs OD', [nexttile(t, 3) nexttile(t, 3+cc)]);
            compute_intersect_gui(onoff_interpolated_smooth, sf_interpolated_smooth, 1, 1, 0, 3, 'jet', 'ONOFF vs SF', [nexttile(t, 9) nexttile(t, 9+cc)]);
            compute_intersect_gui(onoff_interpolated_smooth, ori_map_interpolated, 1, 10, 0, 3, 'jet', 'ONOFF vs Ori', [nexttile(t, 6) nexttile(t, 6+cc)]);

            % sf
            compute_intersect_gui( sf_interpolated_smooth, RetinotopyX_interpolated_smooth, 2, 20, 0, 3, jet, 'SF vs Azimuth', [nexttile(t, 5) nexttile(t, 5+cc)]);
            compute_intersect_gui( sf_interpolated_smooth, RetinotopyY_interpolated_smooth, 2, 20, 0, 3, jet, 'SF vs Elevation', [nexttile(t, 7) nexttile(t, 7+cc)]);

            % Ret
            compute_intersect_gui( RetinotopyX_interpolated_smooth, RetinotopyY_interpolated_smooth, 20, 20, 0, 3, jet, 'Azimuth vs Elevation', [nexttile(t, 2) nexttile(t, 2+cc)]);


        end

        function appdata = updateElectrodArray(app, appdata, def_pos)
            % def_pos   :   default position of electrodes

            n_interpol = appdata.n_interpol;

            appdata = data_based_on_switch(app, app.appdata);

            max_on_response = appdata.plot.max_on_response;
            max_off_response = appdata.plot.max_off_response;
            ODI = appdata.plot.ODI;
            eye_input = appdata.plot.eye_input;

            % plot onoff map
            app.TextArea_processing.Visible = 'off';
            onoff_resp = max_on_response + (-1) * max_off_response;
            onoff_norm = onoff_resp / max(abs(onoff_resp(:)));
            onoff_norm = ODI .* onoff_norm;
            onoff_norm_interp = imresize(onoff_norm, n_interpol);
            axes_onoff = app.UIAxes_rf_select_map;
            axes_onoff.Visible = 'on';
            axes_onoff.Toolbar.Visible = 'off';
            cla(axes_onoff)
            imagesc(axes_onoff, onoff_norm_interp), axis(axes_onoff, 'square'), colormap(axes_onoff,'jet'),
            hold(axes_onoff, 'on')
            contour(axes_onoff, appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', app.od_contour_w);
            axis(axes_onoff, 'off')
            title(axes_onoff, eye_input)
            axes_onoff.XColor = [0 0 0];
            axes_onoff.YColor = [0 0 0];

            app.ElectrodePanel.Title = ['Cortical receptive fields ' '(' eye_input ')'];
            %  title(app.ElectrodePanel.Title, ['Cortical receptive fields' '(' eye_input ')'])

            app.GridLayout13.Visible = 'on';
            app.GridLayout12.Visible = 'on';

            if def_pos < 2
                if def_pos == 1
                    % defualt line
                    map_width = size(app.appdata.ODCrtxPlt_interpolated, 1);
                    map_half_width = round(map_width/2);
                    roi = drawline(axes_onoff, 'Color', 'r', 'Position', [1 map_half_width;map_width map_half_width]);
                else
                    roi = drawline(axes_onoff, 'Color', 'r');
                end

            else
                roi = drawline(axes_onoff, 'Color', 'r', 'Position', [str2double(app.X1_rf.Value) str2double(app.Y1_rf.Value); str2double(app.X2_rf.Value) str2double(app.Y2_rf.Value)]);
            end

            point_pos = roi.Position;
            X1 = point_pos(1, 1);
            Y1 = point_pos(1, 2);
            X2 = point_pos(2, 1);
            Y2 = point_pos(2, 2);
            app.X1_rf.Value = num2str(ceil(X1));
            app.X2_rf.Value = num2str(floor(X2));
            app.Y1_rf.Value = num2str(ceil(Y1));
            app.Y2_rf.Value = num2str(floor(Y2));

            X1 = ceil(str2double(app.X1_rf.Value) / n_interpol);
            X2 = floor(str2double(app.X2_rf.Value) / n_interpol);
            Y1 = ceil(str2double(app.Y1_rf.Value) / n_interpol);
            Y2 = floor(str2double(app.Y2_rf.Value) / n_interpol);
            appdata.plot.ind_row_electrode = ceil(linspace(Y1, Y2, appdata.num_electrode));
            appdata.plot.ind_col_electrode = ceil(linspace(X1, X2, appdata.num_electrode));

            app.appdata = appdata;
            update_rf(app, app.appdata)
        end

        function appdata = data_based_on_switch(app, appdata)
            % determine the data to show based on the eye input switch status
            if app.Switchvalue == 1
                appdata.plot.eye_input = 'Contra';
                appdata.plot.ODI = appdata.data_primord_contra.ODI;

                appdata.plot.max_on_response = appdata.data_contra_mature.max_on_response_norm;
                appdata.plot.max_off_response = appdata.data_contra_mature.max_off_response_norm;

                appdata.plot.allCXrfONReference = appdata.rf_contra_mature.allCXrfON ;
                appdata.plot.allCXrfOFFReference = appdata.rf_contra_mature.allCXrfOFF ;
                appdata.plot.OriHistRFReference =  appdata.data_contra_mature.OriTunPolarPlot;
                appdata.plot.OriBinRFReference =  appdata.data_contra_mature.angles;

                appdata.plot.axis_color = [0 0 0];
            else
                appdata.plot.eye_input = 'Ipsi';
                appdata.plot.ODI = appdata.data_primord_ipsi.ODI ;

                appdata.plot.max_on_response = appdata.data_ipsi_mature.max_on_response_norm;
                appdata.plot.max_off_response = appdata.data_ipsi_mature.max_off_response_norm;

                appdata.plot.allCXrfONReference = appdata.rf_ipsi_mature.allCXrfON ;
                appdata.plot.allCXrfOFFReference = appdata.rf_ipsi_mature.allCXrfOFF ;
                appdata.plot.OriHistRFReference =  appdata.data_ipsi_mature.OriTunPolarPlot;
                appdata.plot.OriBinRFReference =  appdata.data_ipsi_mature.angles;

                appdata.plot.axis_color = [0.9647    0.4235    0.0314]; % orange
            end
            appdata.plot.caxisVal = 1;
            app.appdata = appdata;
        end

        function update_rf(app, appdata)

            ind_row = appdata.plot.ind_row_electrode;
            ind_col = appdata.plot.ind_col_electrode;
            allCXrfONReference = appdata.plot.allCXrfONReference;
            allCXrfOFFReference = appdata.plot.allCXrfOFFReference;
            ODI = appdata.plot.ODI;

            t1 = tiledlayout(app.ElectrodePanel, 2, 10, 'TileSpacing', 'none', 'Padding', 'none');
            % t2 = tiledlayout(app.AfferentconvergenceTab, 2, 10, 'TileSpacing', 'none', 'Padding', 'none');
            
            app.Panel_aff_convergence.Visible = 'on'; 
            t2 = tiledlayout(app.Panel_aff_convergence, 2, 11);

            % t3 = tiledlayout(app.AfferentaxonarborsTab, 2, 10);
            app.Panel_aff_axon_arbors.Visible = 'on'; 
            t3 = tiledlayout(app.Panel_aff_axon_arbors, 2, 11, 'TileSpacing', 'compact', 'Padding', 'loose');

            
            for i = 1:appdata.num_electrode
                ax1 = nexttile(t1);
                colormap(ax1, 'jet');

                RFrefON  =  allCXrfONReference{ind_row(i), ind_col(i)};
                RFrefOFF = allCXrfOFFReference{ind_row(i), ind_col(i)};
                RFreONOFF = RFrefON + RFrefOFF;
                MaxValrf  = max(max(abs(RFrefON(:))), max(abs(RFrefOFF(:))));

                RFreONOFFnorm = RFreONOFF / MaxValrf;
                RFreONOFFnorm = RFreONOFFnorm .* ODI(ind_row(i), ind_col(i));

                plot_rf(app, app.appdata, RFreONOFFnorm, ax1)
                plot_ori_tuning_curve(app, app.appdata, ind_row(i), ind_col(i), nexttile(t1, 10+i))

                ax2 = nexttile(t2);
                ax3 = nexttile(t2, 11+i);

                % colormap(ax2, 'jet');
                % plot_rf(app, app.appdata, RFreONOFFnorm, ax2)


                plot_aff_vis_space_gui(appdata, RFreONOFFnorm, ind_row(i), ind_col(i), ax2, ax3, appdata.plot.axis_color, appdata.plot.eye_input)
                setAxes(app, ax2);
                setAxes(app, ax3);

                % plot TAA (make another function for this part to avoid confution because it is not part of updating rf)
                ax4 = nexttile(t3);
                ax5 = nexttile(t3,  11+i);
                plot_TAA_gui(appdata, ind_row(i), ind_col(i), ax4, ax5, appdata.plot.axis_color, appdata.plot.eye_input)
                setAxes(app, ax4);
                setAxes(app, ax5);
            end
        end

        function plot_rf(app, appdata, RFreONOFFnorm, ax)
            caxisVal = appdata.plot.caxisVal;
            axis_color = appdata.plot.axis_color;
            % cla(ax),
            imagesc(ax, RFreONOFFnorm)
            caxis(ax , [-caxisVal caxisVal])
            axis(ax, 'square')
            colormap(ax, 'jet');
            ax.XColor = axis_color;
            ax.YColor = axis_color;
            setAxes(app, ax);
        end

        function plot_ori_tuning_curve(app, appdata, cortex_row, cortex_col, ax)
            OriBinRFReference = appdata.plot.OriBinRFReference;
            OriHistRFReference = appdata.plot.OriHistRFReference;
            axis_color = appdata.plot.axis_color;
            ODI = appdata.plot.ODI;

            app.ElectrodePanel.Visible = 'on';
            %for i=1:appdata.num_electrode
            % ax = nexttile(t);
            angles     = OriBinRFReference{cortex_row, cortex_col};
            ori_tuning = OriHistRFReference{cortex_row, cortex_col};
            ori_tuning_odi = ori_tuning * ODI(cortex_row, cortex_col);
            polarS3(ax, angles, ori_tuning_odi', axis_color);
            pause(0.001)
            %end
        end

        function updateCIlamps(app)
            if app.Switchvalue == 1
                app.Lamp.Enable = false; % contra
                app.Lamp_2.Enable = true;
            else
                app.Lamp.Enable = true; % ipsi
                app.Lamp_2.Enable = false;
            end
        end

        function clear_figures(app)
            cla(app.UIAxes_RGC_contra)
            cla(app.UIAxes_RGC_ipsi)
            cla(app.UIAxes_aff_RF_vis_space)
            title(app.UIAxes_aff_RF_vis_space, '')

            cla(app.UIAxes_Aff_retinotopy)
            cla(app.UIAxes_Segregation_OD)
            cla(app.UIAxes_OD_ONOFF)

            cla(app.UIAxes_mature_orimap_ipsi)
            cla(app.UIAxes_primord_orimap_ipsi)
            cla(app.UIAxes_mature_orimap_contra)
            cla(app.UIAxes_primord_orimap_contra)
            
            cla(app.UIAxes_ONOFF_contra)
            cla(app.UIAxes_ONOFF_ipsi)
            
%             cla(app.UIAxes3_orimap_mature)
%             cla(app.UIAxes2_ONOFF_mature)
%             cla(app.UIAxes_OD_mature)

            cla(app.UIAxes3_orimap_mature)
            cla(app.UIAxes2_ONOFF_mature)
            cla(app.UIAxes_OD_mature)
            cla(app.UIAxes_retinotopy)

            cla(app.UIAxes_SF50)
            cla(app.UIAxes_LPI)
            cla(app.UIAxes_CV)
            cla(app.UIAxes_LHI)

            cla(app.UIAxes_stat_cv_lhi)
            cla(app.UIAxes_stat_lpi_lhi)
            cla(app.UIAxes_stat_sf_lhi)
            cla(app.UIAxes_stat_lpi_cv)
            cla(app.UIAxes_stat_sf_cv)
            cla(app.UIAxes_slope1)
            cla(app.UIAxes_slope2)

            cla(app.UIAxes_ori_od_polar_hist)
            cla(app.UIAxes_od_biased)
            cla(app.UIAxes_orimap_biased)
            cla(app.UIAxes_ori_od_hist)
            cla(app.UIAxes_orimap_contour)

            cla(app.UIAxes_rf_select_map)

            app.GridLayout_ori_biased.Visible = 'off';
            app.ElectrodePanel.Visible = 'off';
            app.GridLayout12.Visible = 'off';
            app.TextArea_processing.Visible = 'on';
            %app.UIAxes_ASF.Visible = 'off';
            cla(app.UIAxes_ASF)
            %             cla(app.GridLayout_rfs)

            app.Panel_aff_axon_arbors.Visible = 'off'; 
            app.Panel_aff_convergence.Visible = 'off'; 
            app.map_dim_relations2.Visible = 'off'; 
            title(app.UIAxes_rf_select_map, '')
            
            % app.GridLayout_map_dim2.Visible = 'off';
            
            pause(0.01)
        end

        function update_asf_figure(app, appdata)
            ax = app.UIAxes_ASF;
            ax.Visible = 'on';

            cla(ax)
            CenterRadiusX = appdata.CenterRadiusX;
            CenterRadiusY = appdata.CenterRadiusY;
            SurroundRadiusX = appdata.SurroundRadiusX;
            SurroundRadiusY = appdata.SurroundRadiusY;

            CrtxLength = appdata.CrtxLength;
            hh = round(CrtxLength/2); % half cortical length

            %imagesc(ax, zeros(CrtxLength))
            hold(ax,'on')
            rectangle(ax, 'Position', [hh-CenterRadiusX/2, hh-CenterRadiusY/2, CenterRadiusX, CenterRadiusY],'EdgeColor', 'k', 'LineWidth', 2, 'curvature', [1 1])
            rectangle(ax, 'Position', [hh-SurroundRadiusX/2, hh-SurroundRadiusY/2, SurroundRadiusX, SurroundRadiusY],'EdgeColor', 'k', 'LineWidth', 2, 'curvature', [1 1])
            set(ax, 'ydir', 'reverse')

            title(ax, sprintf('Afferent sorting \n filter'))
            axis(ax, 'square')
            %             xlim(ax, [0-10 CrtxLength+10])
            %             ylim(ax, [0-10 CrtxLength+10])
            xlim(ax, [0 CrtxLength])
            ylim(ax, [0 CrtxLength])

            % axis(ax, 'off')
            setAxes(app, ax);
            ax.XColor = [0 0 0];
            ax.YColor = [0 0 0];
            pause(0.001)
        end
    end

    methods (Access = public)
        % Set axis properties
        function setAxes(~,axes)
            axes.Visible = 'on';
            set(axes, 'tickdir', 'out');
            set(axes, 'linewidth', 2);
            set(axes, 'Box', 'on');
            set(axes, 'xtick', []);
            set(axes, 'ytick', []);
            pause(0.001)
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            units = get(groot, 'Units'); % backup Units
            set(groot, 'Units', 'pixels');
            app.UIFigure.Position = get(groot, 'ScreenSize');
            set(groot, 'Units', units); % restore original units

            %             screenSize = get(groot,'ScreenSize');
            %             screenWidth = screenSize(3);
            %             screenHeight = screenSize(4);
            %             left = screenWidth*0.1;
            %             bottom = screenHeight*0.1ElectrodePanel;
            %             width = screenWidth*0.8;
            %             height = screenHeight*0.8;
            %             drawnow;
            %             app.UIFigure.Position = [left bottom width height];
        end

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)

            app.appdata = [];
            addpath([pwd '/functions_user_friendly']);
            clear_figures(app)

            debug = 0;
            show_fig = 0;
            save_file = 0;

            % Loading user's main parameters
            % app.appdata.rng_trial = app.jitter_seed_EditField.Value;
            app.appdata.rng_trial = str2double(app.JitterseedDropDown.Value);

            app.appdata.electrode_position = 10;
            app.appdata.num_electrode = 10;
            app.appdata.alpha = 0.6;     % percent of maximum response for iso-retinotopic measurement
            app.appdata.n_interpol = 4;  % map interpolation for presentation purposes

            app.appdata.od_contour_levels  = 1; %...
            app.appdata.od_contour = 1;

            
            app.appdata.contra_font_color = [0 0 0]/255;
            app.appdata.ipsi_font_color   = [245 134 52]/255;

            % Loading user's retina parameters
            app.appdata.retina.distAPON  = str2double(app.RGCdistxDropDown.Value);
            app.appdata.retina.distMLON  = str2double(app.RGCdistyDropDown.Value);
            app.appdata.retina.distAPOFF = str2double(app.RGCdistxDropDown.Value);
            app.appdata.retina.distMLOFF = str2double(app.RGCdistyDropDown.Value);

            app.appdata.retina.jitterONx  = str2double(app.JitterONxDropDown.Value);
            app.appdata.retina.jitterONy  = str2double(app.JitterONyDropDown.Value);
            app.appdata.retina.jitterOFFx = str2double(app.JitterOFFxDropDown.Value);
            app.appdata.retina.jitterOFFy = str2double(app.JitterOFFyDropDown.Value);

            app.appdata.retina.MinDistONOFF = 1;
            app.appdata.retina.MaxDistONOFF = 3.4;

            % afferent density factor 
            % app.appdata.retina.aff_density_factor_EditField = app.aff_density_factor_EditField.Value;
            ADF = str2double(app.AffdenfactorDropDown.Value);
            app.appdata.retina.aff_density_factor = ADF; 
            app.appdata.pix2mic = 22;

            % Loading user's cortical parameters
            % app.appdata.CrtxLength = app.CortexLengthEditField.Value;
            CortexLengthDropDownValueChanged(app, event)
            app.appdata.CrtxLength = app.appdata.crtx_length;

            %% Retina
            rng(app.appdata.rng_trial)
            app.appdata.CrtxWidth = app.appdata.CrtxLength; % Width of Cortical Plate
            app.appdata.rowRange = 1:app.appdata.CrtxLength;
            app.appdata.colRange = 1:app.appdata.CrtxLength;
            [app.appdata] = make_retina_app3(app.appdata, show_fig);

            % Updating Retina UI
            updatRetinaUI(app, app.appdata, 'contra');
            updatRetinaUI(app, app.appdata, 'ipsi');

            % afferent spread in cortex (there is no aff available with spread 10 for some cortical locations)
            dist_min_RGC = app.appdata.retina.dist_min_RGC ;
            rAffSpread = (10 / dist_min_RGC) * 5 ; % def dist min: 5, def arbor value : 10 
            rAffSpread = rAffSpread * ADF;
            if rAffSpread < 10 
                % if spread is less than 10, there are cases that no afferent spreading at some cortical locations (mostly near the borders)
                rAffSpread =10; 
            end
            app.appdata.rAffSpread = rAffSpread;  % 10;

            app.appdata.pix2mm = 50; %app.PixelsizeEditField.Value;           % each pixel equals 50 microns in cortex
            app.appdata.sigma_LHI_2d = 100; %app.sigma_LHI_2dEditField.Value ;%   in microns
            app.appdata.n_ori_smooth = 4; %app.n_ori_smoothEditField.Value;   %   n orimap interpolation
            % app.appdata.pix2deg = app.PixelindegreeEditField.Value;
            app.appdata.pix2deg = 40;

            % Loading user's afferent sorting parameters
            app.appdata.NSortOD = 10;    % app.ODsortingiterationsEditField.Value;               % Number of iteration for OD sorting
            app.appdata.NSortONOFF = 10; % app.ONOFFsortingiterationsEditField.Value;            % Number of iteration for ON/OFF sorting

            % Sorting filter parameters come from the RGC distribution
            % CenterRadiusX is measured based on RGC distance in make_retina_app function 
            aff_sampling_density        = app.appdata.CenterRadiusX;
            app.appdata.aff_sampling_density = aff_sampling_density; 
            % the following part is moved to make_retina_app function 
%             if aff_sampling_density < 3
%                 % 3 arbitrary number, I could add more options for RGC distance and the
%                 % aff_sampling_density would be lower than 1 automatically,
%                 % but the running time increases for larger initial RGC
%                 % distance because it increases the visual space size
%                 aff_sampling_density = .5 ;
%             end
            
            if aff_sampling_density < 1
                app.appdata.n_interpol = 1; 
                app.appdata.n_ori_smooth = 1;
                app.appdata.od_contour_levels = 0; 
                app.appdata.od_contour = 0; 
            end 

            app.appdata.CenterRadiusX   = aff_sampling_density;
            app.appdata.CenterRadiusY   = aff_sampling_density;
            app.appdata.SurroundRadiusX = app.appdata.SurroundRadiusX;
            app.appdata.SurroundRadiusY = app.appdata.SurroundRadiusY;

            app.appdata.acf_parameters_coverage.acf_center_x = .75 * aff_sampling_density; %15;
            app.appdata.acf_parameters_coverage.acf_center_y = .75 * aff_sampling_density; %15;
            app.appdata.acf_parameters_coverage.acf_surround_x = 1.5 * aff_sampling_density; %30;
            app.appdata.acf_parameters_coverage.acf_surround_y = .75 * aff_sampling_density; %15;
            %   Swindale
            app.appdata.acf_parameters_smooth.acf_center_x =  .75 * aff_sampling_density; %15;
            app.appdata.acf_parameters_smooth.acf_center_y =  .75 * aff_sampling_density; %15;
            app.appdata.acf_parameters_smooth.acf_surround_x = 1.5 * aff_sampling_density; %30;
            app.appdata.acf_parameters_smooth.acf_surround_y = 1.5 * aff_sampling_density; %30;
            app.appdata.CSRatio = 2 ;              % center surround DoG ratio

            %  Other parameters
            app.appdata.OFFONDensityRatio = 1;
            app.appdata.sf_lSamp = 1:50;    % sampling locations in pixel from the center in fft space

            %if strcmp(app.retina_switch.Value, 'symmetrical')
            app.appdata.r_covered_aff = 20; %20; % average radius covered by RFs in cortex in pixel
            %elseif strcmp(app.retina_switch.Value, 'asymmetrical')
            % smaller value will result it more rf with 0 ori
            % preference, (because the are more afferents in horizontal direction than vertical direction)
            %app.appdata.r_covered_aff = 25; %20; % average radius covered by RFs in cortex in pixel
            %end

            % ADF = app.appdata.retina.aff_density_factor_EditField; % afferent density factor
            rf_sd_RGC = ADF * 5; % default val : 5;
            app.appdata.RFONCenter = rf_sd_RGC; % 5;                        % Receptive Field Center
            app.appdata.RFONSurround = round(1.1*app.appdata.RFONCenter);  % Receptive Field Surround
            app.appdata.RFOFFCenter = app.appdata.RFONCenter;
            app.appdata.RFOFFSurround = app.appdata.RFONSurround;
%             NumRandPoint = 200; % to measure correlations
%             boundary = 0.1 ;% remove 10% of the boundary for correlation measurements

            % Making Colormap for retinotopy
            [app.appdata.cmap,app.appdata.RetinotopyRFspace_plot] = make_colormap(app.appdata.rf_space_x, 0);

            % Putting the cells in Matrix  without changing the location (cortical plate)
            [app.appdata] = make_retina_grid_rf_space_app(app.appdata);

            [app.appdata] = make_cortex_app(app.appdata);

            % plot thalamic afferetns in visual space
            update_thalamic_aff_UI(app, app.appdata)

            % OD Sorting
            app.appdata.SortFiltOD = functionDoGfilter(app.appdata.CenterRadiusX, app.appdata.CenterRadiusY, app.appdata.SurroundRadiusX,app.appdata.SurroundRadiusY,app.appdata.CSRatio);
            [app.appdata] = sort_afferent_od2_app(app.appdata);

            % ONOFF Sorting
            app.appdata.SortFiltONOFF = functionDoGfilter(app.appdata.CenterRadiusX, app.appdata.CenterRadiusY, app.appdata.SurroundRadiusY, app.appdata.SurroundRadiusX, app.appdata.CSRatio); % the Y and X is reveresed for ONOFF filter comparing to OD filter
            [app.appdata] = sort_afferent_onoff2_app(app.appdata);

            [app.appdata] = smooth_map_app(app.appdata);
            %             [app.appdata] = smooth_map_app_sn(app, app.appdata);
            [app.appdata] = show_maps_segregation(app, app.appdata);
            update_asf_figure(app, app.appdata)

            % Updating the smoohting map
            % updatSmoothUI(app, app.appdata, app.UIAxes_Vis_space);

            % labeling islands
            [app.appdata.ONOFFODLabelSorted, app.appdata.CenterPinwheelSorted, app.appdata.NumPinwheel, app.appdata.center_island_image] = label_ononff_od_island_gui(app.appdata, debug, 0);

            % synaptic weight based on islands
            app.appdata.weight_map_island = ones(size(app.appdata.ONOFFODLabelSorted)); % it should be tested whether adding weights here helps bring the pw closer to the center of islands
            % weight_map_island = synaptic_weight_island(ONOFFODLabelSorted,NumPinwheel, CenterPinwheelSorted,center_island_image,ODCrtxPlt_smooth,1);

            %   Making the receptive field in Retinal(Visual) Space with Gaussian functions
            app.appdata.RetinaRF = make_rf_retina(app.appdata.TotalRetinaCells, app.appdata.RFONCenter, app.appdata.RFOFFCenter, app.appdata.IndRetinaRfspace, app.appdata.onoff_rf_space);

            % primordial measurements
            [app.appdata] = generate_rf_cortex_primord_app_sn(app.appdata);
            update_primord_orimap_UI(app, app.appdata);

            if aff_sampling_density > 1 
                % modify RF weight based on orimap (pw locations)
                [app.appdata.weight_map_rf_reference_orimap] = synaptic_weight_island2(app.appdata.reference_ori_map, app.appdata.ONOFFODLabelSorted, app.appdata.ODCrtxPlt_smooth, app.appdata.ONOFF_smoothed, 0);
    
                app.appdata.eyepref = 'contra'; % default contra
                [app.appdata] = generate_rf_cortex_mature_app_sn(app.appdata);
    
                app.appdata.eyepref = 'ipsi';   % default contra
                [app.appdata] = generate_rf_cortex_mature_app_sn(app.appdata);
            else 
                app.appdata.rf_contra_mature = app.appdata.rf_primord_contra; 
                app.appdata.data_contra_mature = app.appdata.data_primord_contra; 
                
                app.appdata.rf_ipsi_mature = app.appdata.rf_primord_ipsi; 
                app.appdata.data_ipsi_mature = app.appdata.data_primord_ipsi; 
                
                [app.appdata] = primord_to_mature_data_transfer_gui(app.appdata);
            end 
            
            update_mature_orimaps_UI(app, app.appdata);

            update_mature_onoff_maps_UI(app, app.appdata, 'contra')
            update_mature_onoff_maps_UI(app, app.appdata, 'ipsi')

            % analysis_pref_eye_uf('contra', data_contra_mature, rf_contra_mature, ODCrtxPlt_smooth, ODCrtxPlt_interpolated, electrode_position, boundary, NumRandPoint, 1)
            update_mature_maps_UI(app, app.appdata);

            % app.GridLayout_map_dim2.Visible = 'on';
            update_maps_intersection_UI(app, app.appdata)

            generate_fig4_contra_ipsi_app(app, app.appdata);

            % generate the orientation map, ocular dominance contours' intersection
            analysis_ori_map_app(app, app.appdata)

            % update array presentation
            updateElectrodArray(app, app.appdata, 1);

            app.GridLayout_ori_biased.Visible = 'on';
        end

        % Value changed function: EyeofpreferenceSwitch
        function EyeofpreferenceSwitchValueChanged(app, event)
            if strcmp(app.EyeofpreferenceSwitch.Value, 'Contra')
                app.Switchvalue = 1;
            else
                app.Switchvalue = 2;
            end
            updateCIlamps(app);
            updateElectrodArray(app, app.appdata, 2);
        end

        % Window button down function: UIFigure
        function UIFigureWindowButtonDown(app, event)
            %             temp = app.UIAxes_Segregation_OD.CurrentPoint; % Returns 2x3 array of points
            %             loc = [temp(1,1) temp(1,2)]; % Gets the (x,y) coordinates
            %
            %             if app.CurrentPoint ==1
            %                 app.FirstPoint = ceil(loc);
            %                 app.CurrentPoint = 2;
            %             elseif app.CurrentPoint ==2
            %                 app.SecondPoint = ceil(loc);
            %                 app.CurrentPoint = 1;
            %             end
            %
            %             if app.FinalPrecessDone
            %                 updateElectrodArray(app,app.appdata);
            %             end
        end

        % Callback function: CloseButton, UIAxes_RGC_contra
        function CloseButtonPushed(app, event)
            closereq
        end

        % Button pushed function: DrawalineonthemapButton
        function DrawalineonthemapButtonPushed(app, event)
            %             ax = app.UIAxes_rf_select_map;
            %             map_width = size(app.appdata.ODCrtxPlt_interpolated,1);
            %             map_half_width = round(map_width/2);
            %             roi = drawline(ax, 'Color', 'r', 'Position', [1 map_half_width;map_width map_half_width]);
            %             point_pos = roi.Position;
            %             X1 = point_pos(1, 1);
            %             Y1 = point_pos(1, 2);
            %             X2 = point_pos(2, 1);
            %             Y2 = point_pos(2, 2);
            %             app.X1_rf.Value = num2str(ceil(X1)) ;
            %             app.X2_rf.Value = num2str(floor(X2)) ;
            %             app.Y1_rf.Value = num2str(ceil(Y1)) ;
            %             app.Y2_rf.Value = num2str(floor(Y2)) ;

            updateElectrodArray(app, app.appdata, 0);




        end

        % Button pushed function: GenerateButton
        function GenerateButtonPushed(app, event)
            app.appdata.ori_biased = app.OrientationbiasedEditField.Value;
            %app.appdata.ori_biased_amount = app.BiasamountEditField.Value;
            app.appdata.ori_biased_amount = str2double(app.BiaspercentDropDown.Value)/100;

            %app.UIAxes_orimap_biased.Visible = 'off';
            cla(app.UIAxes_orimap_biased)
            app.TextArea_processing.Visible = 'on';
            pause(0.001)
            generate_rf_cortex_ori_biased_app(app, app.appdata)
            app.TextArea_processing.Visible = 'off';
        end

        % Value changed function: CortexLengthDropDown
        function CortexLengthDropDownValueChanged(app, event)
            value = app.CortexLengthDropDown.Value;
            if strcmp(value, 'Cortical area: 2*2 mm')
                crtx_length = 40;
            elseif strcmp(value, 'Cortical area: 3*3 mm')
                crtx_length = 60;
            end
            app.appdata.crtx_length = crtx_length;
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1412 699];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.Resize = 'off';
            app.UIFigure.WindowButtonDownFcn = createCallbackFcn(app, @UIFigureWindowButtonDown, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {20, 380, 10, '1x', 10, '1x', 10, '1x', 10, '1x', 10, '1x', 10};
            app.GridLayout.RowHeight = {'1x', 20, '1x', '1x', 30};

            % Create TabGroup
            app.TabGroup = uitabgroup(app.GridLayout);
            app.TabGroup.Layout.Row = 1;
            app.TabGroup.Layout.Column = 2;

            % Create MainTab
            app.MainTab = uitab(app.TabGroup);
            app.MainTab.Title = 'Main';

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.MainTab);
            app.GridLayout3.ColumnWidth = {'1x'};
            app.GridLayout3.RowHeight = {'1x', '1x', '1x', '1x'};

            % Create StartButton
            app.StartButton = uibutton(app.GridLayout3, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton.FontWeight = 'bold';
            app.StartButton.Layout.Row = 2;
            app.StartButton.Layout.Column = 1;
            app.StartButton.Text = 'Start';

            % Create CloseButton
            app.CloseButton = uibutton(app.GridLayout3, 'push');
            app.CloseButton.ButtonPushedFcn = createCallbackFcn(app, @CloseButtonPushed, true);
            app.CloseButton.Layout.Row = 3;
            app.CloseButton.Layout.Column = 1;
            app.CloseButton.Text = 'Close ';

            % Create CortexLengthDropDown
            app.CortexLengthDropDown = uidropdown(app.GridLayout3);
            app.CortexLengthDropDown.Items = {'Cortical area: 2*2 mm', 'Cortical area: 3*3 mm', ''};
            app.CortexLengthDropDown.ValueChangedFcn = createCallbackFcn(app, @CortexLengthDropDownValueChanged, true);
            app.CortexLengthDropDown.Layout.Row = 1;
            app.CortexLengthDropDown.Layout.Column = 1;
            app.CortexLengthDropDown.Value = 'Cortical area: 2*2 mm';

            % Create TextArea_processing
            app.TextArea_processing = uitextarea(app.GridLayout3);
            app.TextArea_processing.Editable = 'off';
            app.TextArea_processing.HorizontalAlignment = 'center';
            app.TextArea_processing.FontSize = 16;
            app.TextArea_processing.FontColor = [1 0 0];
            app.TextArea_processing.Visible = 'off';
            app.TextArea_processing.Layout.Row = 4;
            app.TextArea_processing.Layout.Column = 1;
            app.TextArea_processing.Value = {'Processing ...'};

            % Create RetinathalamusTab
            app.RetinathalamusTab = uitab(app.TabGroup);
            app.RetinathalamusTab.Title = 'Retina-thalamus';

            % Create GridLayout20
            app.GridLayout20 = uigridlayout(app.RetinathalamusTab);
            app.GridLayout20.ColumnWidth = {'1x', '0.75x', '0.1x', '1x', '0.75x'};
            app.GridLayout20.RowHeight = {'1x', '1x', '1x', '1x'};

            % Create JitterONxDropDownLabel
            app.JitterONxDropDownLabel = uilabel(app.GridLayout20);
            app.JitterONxDropDownLabel.Layout.Row = 1;
            app.JitterONxDropDownLabel.Layout.Column = 4;
            app.JitterONxDropDownLabel.Text = 'Jitter ON x';

            % Create JitterONxDropDown
            app.JitterONxDropDown = uidropdown(app.GridLayout20);
            app.JitterONxDropDown.Items = {'0.05', '0.1', '0.15', '0.2', '0.3', '0.4', '0.5', '1'};
            app.JitterONxDropDown.Layout.Row = 1;
            app.JitterONxDropDown.Layout.Column = 5;
            app.JitterONxDropDown.Value = '0.1';

            % Create JitterONyDropDownLabel
            app.JitterONyDropDownLabel = uilabel(app.GridLayout20);
            app.JitterONyDropDownLabel.Layout.Row = 2;
            app.JitterONyDropDownLabel.Layout.Column = 4;
            app.JitterONyDropDownLabel.Text = 'Jitter ON y';

            % Create JitterONyDropDown
            app.JitterONyDropDown = uidropdown(app.GridLayout20);
            app.JitterONyDropDown.Items = {'0.05', '0.1', '0.15', '0.2', '0.3', '0.4', '0.5', '1'};
            app.JitterONyDropDown.Layout.Row = 2;
            app.JitterONyDropDown.Layout.Column = 5;
            app.JitterONyDropDown.Value = '0.1';

            % Create JitterOFFxDropDownLabel
            app.JitterOFFxDropDownLabel = uilabel(app.GridLayout20);
            app.JitterOFFxDropDownLabel.Layout.Row = 3;
            app.JitterOFFxDropDownLabel.Layout.Column = 4;
            app.JitterOFFxDropDownLabel.Text = 'Jitter OFF x';

            % Create JitterOFFxDropDown
            app.JitterOFFxDropDown = uidropdown(app.GridLayout20);
            app.JitterOFFxDropDown.Items = {'0.05', '0.1', '0.15', '0.2', '0.3', '0.4', '0.5', '1'};
            app.JitterOFFxDropDown.Layout.Row = 3;
            app.JitterOFFxDropDown.Layout.Column = 5;
            app.JitterOFFxDropDown.Value = '0.15';

            % Create JitterOFFyDropDownLabel
            app.JitterOFFyDropDownLabel = uilabel(app.GridLayout20);
            app.JitterOFFyDropDownLabel.Layout.Row = 4;
            app.JitterOFFyDropDownLabel.Layout.Column = 4;
            app.JitterOFFyDropDownLabel.Text = 'Jitter OFF y';

            % Create JitterOFFyDropDown
            app.JitterOFFyDropDown = uidropdown(app.GridLayout20);
            app.JitterOFFyDropDown.Items = {'0.05', '0.1', '0.15', '0.2', '0.3', '0.4', '0.5', '1'};
            app.JitterOFFyDropDown.Layout.Row = 4;
            app.JitterOFFyDropDown.Layout.Column = 5;
            app.JitterOFFyDropDown.Value = '0.15';

            % Create AffdenfactorDropDownLabel
            app.AffdenfactorDropDownLabel = uilabel(app.GridLayout20);
            app.AffdenfactorDropDownLabel.Layout.Row = 2;
            app.AffdenfactorDropDownLabel.Layout.Column = 1;
            app.AffdenfactorDropDownLabel.Text = 'Aff den factor';

            % Create AffdenfactorDropDown
            app.AffdenfactorDropDown = uidropdown(app.GridLayout20);
            app.AffdenfactorDropDown.Items = {'1', '2'};
            app.AffdenfactorDropDown.Layout.Row = 2;
            app.AffdenfactorDropDown.Layout.Column = 2;
            app.AffdenfactorDropDown.Value = '1';

            % Create RGCdistxDropDownLabel
            app.RGCdistxDropDownLabel = uilabel(app.GridLayout20);
            app.RGCdistxDropDownLabel.Layout.Row = 3;
            app.RGCdistxDropDownLabel.Layout.Column = 1;
            app.RGCdistxDropDownLabel.Text = 'RGC dist x';

            % Create RGCdistxDropDown
            app.RGCdistxDropDown = uidropdown(app.GridLayout20);
            app.RGCdistxDropDown.Items = {'3', '4', '5'};
            app.RGCdistxDropDown.Layout.Row = 3;
            app.RGCdistxDropDown.Layout.Column = 2;
            app.RGCdistxDropDown.Value = '5';

            % Create RGCdistyDropDownLabel
            app.RGCdistyDropDownLabel = uilabel(app.GridLayout20);
            app.RGCdistyDropDownLabel.Layout.Row = 4;
            app.RGCdistyDropDownLabel.Layout.Column = 1;
            app.RGCdistyDropDownLabel.Text = 'RGC dist y';

            % Create RGCdistyDropDown
            app.RGCdistyDropDown = uidropdown(app.GridLayout20);
            app.RGCdistyDropDown.Items = {'3', '4', '5'};
            app.RGCdistyDropDown.Layout.Row = 4;
            app.RGCdistyDropDown.Layout.Column = 2;
            app.RGCdistyDropDown.Value = '5';

            % Create JitterseedDropDownLabel
            app.JitterseedDropDownLabel = uilabel(app.GridLayout20);
            app.JitterseedDropDownLabel.Layout.Row = 1;
            app.JitterseedDropDownLabel.Layout.Column = 1;
            app.JitterseedDropDownLabel.Text = 'Jitter seed';

            % Create JitterseedDropDown
            app.JitterseedDropDown = uidropdown(app.GridLayout20);
            app.JitterseedDropDown.Items = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
            app.JitterseedDropDown.Layout.Row = 1;
            app.JitterseedDropDown.Layout.Column = 2;
            app.JitterseedDropDown.Value = '1';

            % Create GridLayout8
            app.GridLayout8 = uigridlayout(app.GridLayout);
            app.GridLayout8.RowHeight = {'1x', '1x', '1x', '1x'};
            app.GridLayout8.Layout.Row = 3;
            app.GridLayout8.Layout.Column = 2;

            % Create ElectrodePanel
            app.ElectrodePanel = uipanel(app.GridLayout);
            app.ElectrodePanel.Title = 'Cortical receptive fields';
            app.ElectrodePanel.Layout.Row = 4;
            app.ElectrodePanel.Layout.Column = [3 12];

            % Create GridLayout11
            app.GridLayout11 = uigridlayout(app.GridLayout);
            app.GridLayout11.RowHeight = {'1x'};
            app.GridLayout11.Layout.Row = 4;
            app.GridLayout11.Layout.Column = 2;

            % Create UIAxes_rf_select_map
            app.UIAxes_rf_select_map = uiaxes(app.GridLayout11);
            zlabel(app.UIAxes_rf_select_map, 'Z')
            app.UIAxes_rf_select_map.Toolbar.Visible = 'off';
            app.UIAxes_rf_select_map.Layout.Row = 1;
            app.UIAxes_rf_select_map.Layout.Column = 1;
            app.UIAxes_rf_select_map.Visible = 'off';

            % Create GridLayout12
            app.GridLayout12 = uigridlayout(app.GridLayout11);
            app.GridLayout12.ColumnWidth = {'1x'};
            app.GridLayout12.Layout.Row = 1;
            app.GridLayout12.Layout.Column = 2;

            % Create GridLayout13
            app.GridLayout13 = uigridlayout(app.GridLayout12);
            app.GridLayout13.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout13.RowHeight = {'1x', '1x', '1x', '1x'};
            app.GridLayout13.Visible = 'off';
            app.GridLayout13.Layout.Row = [1 2];
            app.GridLayout13.Layout.Column = 1;

            % Create DrawalineonthemapButton
            app.DrawalineonthemapButton = uibutton(app.GridLayout13, 'push');
            app.DrawalineonthemapButton.ButtonPushedFcn = createCallbackFcn(app, @DrawalineonthemapButtonPushed, true);
            app.DrawalineonthemapButton.FontSize = 11;
            app.DrawalineonthemapButton.Layout.Row = 1;
            app.DrawalineonthemapButton.Layout.Column = [1 6];
            app.DrawalineonthemapButton.Text = 'Draw a line on the map';

            % Create TextArea
            app.TextArea = uitextarea(app.GridLayout13);
            app.TextArea.Editable = 'off';
            app.TextArea.FontSize = 11;
            app.TextArea.BackgroundColor = [0.902 0.902 0.902];
            app.TextArea.Visible = 'off';
            app.TextArea.Layout.Row = 2;
            app.TextArea.Layout.Column = [1 6];
            app.TextArea.Value = {'Draw a line on the Map'};

            % Create EyeofpreferenceSwitch
            app.EyeofpreferenceSwitch = uiswitch(app.GridLayout13, 'toggle');
            app.EyeofpreferenceSwitch.Items = {'Contra', 'Ipsi'};
            app.EyeofpreferenceSwitch.Orientation = 'horizontal';
            app.EyeofpreferenceSwitch.ValueChangedFcn = createCallbackFcn(app, @EyeofpreferenceSwitchValueChanged, true);
            app.EyeofpreferenceSwitch.FontSize = 11;
            app.EyeofpreferenceSwitch.Layout.Row = 4;
            app.EyeofpreferenceSwitch.Layout.Column = [1 6];
            app.EyeofpreferenceSwitch.Value = 'Contra';

            % Create EyeofpreferenceSwitchLabel
            app.EyeofpreferenceSwitchLabel = uilabel(app.GridLayout13);
            app.EyeofpreferenceSwitchLabel.HorizontalAlignment = 'center';
            app.EyeofpreferenceSwitchLabel.FontSize = 11;
            app.EyeofpreferenceSwitchLabel.Layout.Row = 3;
            app.EyeofpreferenceSwitchLabel.Layout.Column = [2 5];
            app.EyeofpreferenceSwitchLabel.Text = 'Eye of preference';

            % Create Lamp_2
            app.Lamp_2 = uilamp(app.GridLayout13);
            app.Lamp_2.Layout.Row = 3;
            app.Lamp_2.Layout.Column = 1;

            % Create Lamp
            app.Lamp = uilamp(app.GridLayout13);
            app.Lamp.Enable = 'off';
            app.Lamp.Layout.Row = 3;
            app.Lamp.Layout.Column = 6;

            % Create GridLayout14
            app.GridLayout14 = uigridlayout(app.GridLayout);
            app.GridLayout14.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout14.RowHeight = {'1x'};
            app.GridLayout14.Layout.Row = 5;
            app.GridLayout14.Layout.Column = 2;

            % Create X1_rf
            app.X1_rf = uieditfield(app.GridLayout14, 'text');
            app.X1_rf.Editable = 'off';
            app.X1_rf.Visible = 'off';
            app.X1_rf.Layout.Row = 1;
            app.X1_rf.Layout.Column = 2;

            % Create X1EditFieldLabel
            app.X1EditFieldLabel = uilabel(app.GridLayout14);
            app.X1EditFieldLabel.HorizontalAlignment = 'right';
            app.X1EditFieldLabel.Visible = 'off';
            app.X1EditFieldLabel.Layout.Row = 1;
            app.X1EditFieldLabel.Layout.Column = 1;
            app.X1EditFieldLabel.Text = 'X1';

            % Create X2_rf
            app.X2_rf = uieditfield(app.GridLayout14, 'text');
            app.X2_rf.Editable = 'off';
            app.X2_rf.Visible = 'off';
            app.X2_rf.Layout.Row = 1;
            app.X2_rf.Layout.Column = 6;

            % Create X2EditFieldLabel
            app.X2EditFieldLabel = uilabel(app.GridLayout14);
            app.X2EditFieldLabel.HorizontalAlignment = 'right';
            app.X2EditFieldLabel.Visible = 'off';
            app.X2EditFieldLabel.Layout.Row = 1;
            app.X2EditFieldLabel.Layout.Column = 5;
            app.X2EditFieldLabel.Text = 'X2';

            % Create Y1_rf
            app.Y1_rf = uieditfield(app.GridLayout14, 'text');
            app.Y1_rf.Editable = 'off';
            app.Y1_rf.Visible = 'off';
            app.Y1_rf.Layout.Row = 1;
            app.Y1_rf.Layout.Column = 4;

            % Create Y1EditFieldLabel
            app.Y1EditFieldLabel = uilabel(app.GridLayout14);
            app.Y1EditFieldLabel.HorizontalAlignment = 'right';
            app.Y1EditFieldLabel.Visible = 'off';
            app.Y1EditFieldLabel.Layout.Row = 1;
            app.Y1EditFieldLabel.Layout.Column = 3;
            app.Y1EditFieldLabel.Text = 'Y1';

            % Create Y2_rf
            app.Y2_rf = uieditfield(app.GridLayout14, 'text');
            app.Y2_rf.Editable = 'off';
            app.Y2_rf.Visible = 'off';
            app.Y2_rf.Layout.Row = 1;
            app.Y2_rf.Layout.Column = 8;

            % Create Y2EditFieldLabel
            app.Y2EditFieldLabel = uilabel(app.GridLayout14);
            app.Y2EditFieldLabel.HorizontalAlignment = 'right';
            app.Y2EditFieldLabel.Visible = 'off';
            app.Y2EditFieldLabel.Layout.Row = 1;
            app.Y2EditFieldLabel.Layout.Column = 7;
            app.Y2EditFieldLabel.Text = 'Y2';

            % Create TabGroup2
            app.TabGroup2 = uitabgroup(app.GridLayout);
            app.TabGroup2.Layout.Row = [2 3];
            app.TabGroup2.Layout.Column = [2 12];

            % Create CorticaldevelopmentTab
            app.CorticaldevelopmentTab = uitab(app.TabGroup2);
            app.CorticaldevelopmentTab.Title = 'Cortical development';

            % Create GridLayout6
            app.GridLayout6 = uigridlayout(app.CorticaldevelopmentTab);
            app.GridLayout6.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout6.RowHeight = {'0.2x', '1x'};

            % Create UIAxes_mature_orimap_contra
            app.UIAxes_mature_orimap_contra = uiaxes(app.GridLayout6);
            app.UIAxes_mature_orimap_contra.Toolbar.Visible = 'off';
            app.UIAxes_mature_orimap_contra.Layout.Row = 2;
            app.UIAxes_mature_orimap_contra.Layout.Column = 8;
            app.UIAxes_mature_orimap_contra.Visible = 'off';

            % Create UIAxes_primord_orimap_ipsi
            app.UIAxes_primord_orimap_ipsi = uiaxes(app.GridLayout6);
            app.UIAxes_primord_orimap_ipsi.Toolbar.Visible = 'off';
            app.UIAxes_primord_orimap_ipsi.Layout.Row = 2;
            app.UIAxes_primord_orimap_ipsi.Layout.Column = 7;
            app.UIAxes_primord_orimap_ipsi.Visible = 'off';

            % Create UIAxes_mature_orimap_ipsi
            app.UIAxes_mature_orimap_ipsi = uiaxes(app.GridLayout6);
            app.UIAxes_mature_orimap_ipsi.Toolbar.Visible = 'off';
            app.UIAxes_mature_orimap_ipsi.Layout.Row = 2;
            app.UIAxes_mature_orimap_ipsi.Layout.Column = 9;
            app.UIAxes_mature_orimap_ipsi.Visible = 'off';

            % Create UIAxes_Segregation_OD
            app.UIAxes_Segregation_OD = uiaxes(app.GridLayout6);
            app.UIAxes_Segregation_OD.Toolbar.Visible = 'off';
            app.UIAxes_Segregation_OD.Layout.Row = 2;
            app.UIAxes_Segregation_OD.Layout.Column = 1;
            app.UIAxes_Segregation_OD.Visible = 'off';

            % Create UIAxes_Aff_retinotopy
            app.UIAxes_Aff_retinotopy = uiaxes(app.GridLayout6);
            zlabel(app.UIAxes_Aff_retinotopy, 'Z')
            app.UIAxes_Aff_retinotopy.Toolbar.Visible = 'off';
            app.UIAxes_Aff_retinotopy.Layout.Row = 2;
            app.UIAxes_Aff_retinotopy.Layout.Column = 3;
            app.UIAxes_Aff_retinotopy.Visible = 'off';

            % Create UIAxes_OD_ONOFF
            app.UIAxes_OD_ONOFF = uiaxes(app.GridLayout6);
            app.UIAxes_OD_ONOFF.Toolbar.Visible = 'off';
            app.UIAxes_OD_ONOFF.XTickLabelRotation = 0;
            app.UIAxes_OD_ONOFF.YTickLabelRotation = 0;
            app.UIAxes_OD_ONOFF.ZTickLabelRotation = 0;
            app.UIAxes_OD_ONOFF.Layout.Row = 2;
            app.UIAxes_OD_ONOFF.Layout.Column = 2;
            app.UIAxes_OD_ONOFF.Visible = 'off';

            % Create UIAxes_ONOFF_contra
            app.UIAxes_ONOFF_contra = uiaxes(app.GridLayout6);
            app.UIAxes_ONOFF_contra.Toolbar.Visible = 'off';
            app.UIAxes_ONOFF_contra.Layout.Row = 2;
            app.UIAxes_ONOFF_contra.Layout.Column = 4;
            app.UIAxes_ONOFF_contra.Visible = 'off';

            % Create UIAxes_ONOFF_ipsi
            app.UIAxes_ONOFF_ipsi = uiaxes(app.GridLayout6);
            app.UIAxes_ONOFF_ipsi.Toolbar.Visible = 'off';
            app.UIAxes_ONOFF_ipsi.Layout.Row = 2;
            app.UIAxes_ONOFF_ipsi.Layout.Column = 5;
            app.UIAxes_ONOFF_ipsi.Visible = 'off';

            % Create UIAxes_primord_orimap_contra
            app.UIAxes_primord_orimap_contra = uiaxes(app.GridLayout6);
            app.UIAxes_primord_orimap_contra.Toolbar.Visible = 'off';
            app.UIAxes_primord_orimap_contra.Layout.Row = 2;
            app.UIAxes_primord_orimap_contra.Layout.Column = 6;
            app.UIAxes_primord_orimap_contra.Visible = 'off';

            % Create title_eye_pol
            app.title_eye_pol = uieditfield(app.GridLayout6, 'text');
            app.title_eye_pol.Editable = 'off';
            app.title_eye_pol.HorizontalAlignment = 'center';
            app.title_eye_pol.FontSize = 10;
            app.title_eye_pol.Visible = 'off';
            app.title_eye_pol.Layout.Row = 1;
            app.title_eye_pol.Layout.Column = [4 5];
            app.title_eye_pol.Value = 'Eye-polarity grid';

            % Create title_primord_orimap
            app.title_primord_orimap = uieditfield(app.GridLayout6, 'text');
            app.title_primord_orimap.Editable = 'off';
            app.title_primord_orimap.HorizontalAlignment = 'center';
            app.title_primord_orimap.FontSize = 10;
            app.title_primord_orimap.Visible = 'off';
            app.title_primord_orimap.Layout.Row = 1;
            app.title_primord_orimap.Layout.Column = [6 7];
            app.title_primord_orimap.Value = 'Primordial orientation map';

            % Create title_mature_orimap
            app.title_mature_orimap = uieditfield(app.GridLayout6, 'text');
            app.title_mature_orimap.Editable = 'off';
            app.title_mature_orimap.HorizontalAlignment = 'center';
            app.title_mature_orimap.FontSize = 10;
            app.title_mature_orimap.Visible = 'off';
            app.title_mature_orimap.Layout.Row = 1;
            app.title_mature_orimap.Layout.Column = [8 9];
            app.title_mature_orimap.Value = 'Adult orientation map';

            % Create title_segregation_od
            app.title_segregation_od = uieditfield(app.GridLayout6, 'text');
            app.title_segregation_od.Editable = 'off';
            app.title_segregation_od.HorizontalAlignment = 'center';
            app.title_segregation_od.FontSize = 10;
            app.title_segregation_od.Visible = 'off';
            app.title_segregation_od.Layout.Row = 1;
            app.title_segregation_od.Layout.Column = 1;
            app.title_segregation_od.Value = 'Ocular dominance (OD)';

            % Create title_eye_polarity
            app.title_eye_polarity = uieditfield(app.GridLayout6, 'text');
            app.title_eye_polarity.Editable = 'off';
            app.title_eye_polarity.HorizontalAlignment = 'center';
            app.title_eye_polarity.FontSize = 10;
            app.title_eye_polarity.Visible = 'off';
            app.title_eye_polarity.Layout.Row = 1;
            app.title_eye_polarity.Layout.Column = 2;
            app.title_eye_polarity.Value = 'Eye-polarity';

            % Create title_retinotopy
            app.title_retinotopy = uieditfield(app.GridLayout6, 'text');
            app.title_retinotopy.Editable = 'off';
            app.title_retinotopy.HorizontalAlignment = 'center';
            app.title_retinotopy.FontSize = 10;
            app.title_retinotopy.Visible = 'off';
            app.title_retinotopy.Layout.Row = 1;
            app.title_retinotopy.Layout.Column = 3;
            app.title_retinotopy.Value = 'Retinotopy';

            % Create MaturevisualcortexTab
            app.MaturevisualcortexTab = uitab(app.TabGroup2);
            app.MaturevisualcortexTab.Title = 'Mature visual cortex ';

            % Create GridLayout10
            app.GridLayout10 = uigridlayout(app.MaturevisualcortexTab);
            app.GridLayout10.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout10.RowHeight = {'1x'};

            % Create UIAxes_OD_mature
            app.UIAxes_OD_mature = uiaxes(app.GridLayout10);
            app.UIAxes_OD_mature.Toolbar.Visible = 'off';
            app.UIAxes_OD_mature.XTickLabelRotation = 0;
            app.UIAxes_OD_mature.YTickLabelRotation = 0;
            app.UIAxes_OD_mature.ZTickLabelRotation = 0;
            app.UIAxes_OD_mature.Layout.Row = 1;
            app.UIAxes_OD_mature.Layout.Column = 2;
            app.UIAxes_OD_mature.Visible = 'off';

            % Create UIAxes2_ONOFF_mature
            app.UIAxes2_ONOFF_mature = uiaxes(app.GridLayout10);
            app.UIAxes2_ONOFF_mature.Toolbar.Visible = 'off';
            app.UIAxes2_ONOFF_mature.XTickLabelRotation = 0;
            app.UIAxes2_ONOFF_mature.YTickLabelRotation = 0;
            app.UIAxes2_ONOFF_mature.ZTickLabelRotation = 0;
            app.UIAxes2_ONOFF_mature.Layout.Row = 1;
            app.UIAxes2_ONOFF_mature.Layout.Column = 3;
            app.UIAxes2_ONOFF_mature.Visible = 'off';

            % Create UIAxes3_orimap_mature
            app.UIAxes3_orimap_mature = uiaxes(app.GridLayout10);
            app.UIAxes3_orimap_mature.Toolbar.Visible = 'off';
            app.UIAxes3_orimap_mature.XTickLabelRotation = 0;
            app.UIAxes3_orimap_mature.YTickLabelRotation = 0;
            app.UIAxes3_orimap_mature.ZTickLabelRotation = 0;
            app.UIAxes3_orimap_mature.Layout.Row = 1;
            app.UIAxes3_orimap_mature.Layout.Column = 4;
            app.UIAxes3_orimap_mature.Visible = 'off';

            % Create UIAxes_CV
            app.UIAxes_CV = uiaxes(app.GridLayout10);
            app.UIAxes_CV.Toolbar.Visible = 'off';
            app.UIAxes_CV.Layout.Row = 1;
            app.UIAxes_CV.Layout.Column = 5;
            app.UIAxes_CV.Visible = 'off';

            % Create UIAxes_SF50
            app.UIAxes_SF50 = uiaxes(app.GridLayout10);
            app.UIAxes_SF50.Toolbar.Visible = 'off';
            app.UIAxes_SF50.Layout.Row = 1;
            app.UIAxes_SF50.Layout.Column = 7;
            app.UIAxes_SF50.Visible = 'off';

            % Create UIAxes_retinotopy
            app.UIAxes_retinotopy = uiaxes(app.GridLayout10);
            app.UIAxes_retinotopy.Toolbar.Visible = 'off';
            app.UIAxes_retinotopy.Layout.Row = 1;
            app.UIAxes_retinotopy.Layout.Column = 1;
            app.UIAxes_retinotopy.Visible = 'off';

            % Create UIAxes_LHI
            app.UIAxes_LHI = uiaxes(app.GridLayout10);
            app.UIAxes_LHI.Toolbar.Visible = 'off';
            app.UIAxes_LHI.Layout.Row = 1;
            app.UIAxes_LHI.Layout.Column = 6;
            app.UIAxes_LHI.Visible = 'off';

            % Create UIAxes_LPI
            app.UIAxes_LPI = uiaxes(app.GridLayout10);
            app.UIAxes_LPI.Toolbar.Visible = 'off';
            app.UIAxes_LPI.Layout.Row = 1;
            app.UIAxes_LPI.Layout.Column = 8;
            app.UIAxes_LPI.Visible = 'off';

            % Create Mapdimensionrelations1Tab
            app.Mapdimensionrelations1Tab = uitab(app.TabGroup2);
            app.Mapdimensionrelations1Tab.Title = 'Map-dimension relations 1';

            % Create GridLayout15
            app.GridLayout15 = uigridlayout(app.Mapdimensionrelations1Tab);
            app.GridLayout15.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout15.RowHeight = {'1x'};

            % Create UIAxes_stat_sf_cv
            app.UIAxes_stat_sf_cv = uiaxes(app.GridLayout15);
            app.UIAxes_stat_sf_cv.Toolbar.Visible = 'off';
            app.UIAxes_stat_sf_cv.Layout.Row = 1;
            app.UIAxes_stat_sf_cv.Layout.Column = 1;
            app.UIAxes_stat_sf_cv.Visible = 'off';

            % Create UIAxes_stat_lpi_cv
            app.UIAxes_stat_lpi_cv = uiaxes(app.GridLayout15);
            app.UIAxes_stat_lpi_cv.Toolbar.Visible = 'off';
            app.UIAxes_stat_lpi_cv.Layout.Row = 1;
            app.UIAxes_stat_lpi_cv.Layout.Column = 2;
            app.UIAxes_stat_lpi_cv.Visible = 'off';

            % Create UIAxes_stat_sf_lhi
            app.UIAxes_stat_sf_lhi = uiaxes(app.GridLayout15);
            app.UIAxes_stat_sf_lhi.Toolbar.Visible = 'off';
            app.UIAxes_stat_sf_lhi.Layout.Row = 1;
            app.UIAxes_stat_sf_lhi.Layout.Column = 3;
            app.UIAxes_stat_sf_lhi.Visible = 'off';

            % Create UIAxes_stat_lpi_lhi
            app.UIAxes_stat_lpi_lhi = uiaxes(app.GridLayout15);
            app.UIAxes_stat_lpi_lhi.Toolbar.Visible = 'off';
            app.UIAxes_stat_lpi_lhi.Layout.Row = 1;
            app.UIAxes_stat_lpi_lhi.Layout.Column = 4;
            app.UIAxes_stat_lpi_lhi.Visible = 'off';

            % Create UIAxes_stat_cv_lhi
            app.UIAxes_stat_cv_lhi = uiaxes(app.GridLayout15);
            app.UIAxes_stat_cv_lhi.Toolbar.Visible = 'off';
            app.UIAxes_stat_cv_lhi.Layout.Row = 1;
            app.UIAxes_stat_cv_lhi.Layout.Column = 5;
            app.UIAxes_stat_cv_lhi.Visible = 'off';

            % Create UIAxes_slope1
            app.UIAxes_slope1 = uiaxes(app.GridLayout15);
            app.UIAxes_slope1.Toolbar.Visible = 'off';
            app.UIAxes_slope1.Layout.Row = 1;
            app.UIAxes_slope1.Layout.Column = [6 7];
            app.UIAxes_slope1.Visible = 'off';

            % Create UIAxes_slope2
            app.UIAxes_slope2 = uiaxes(app.GridLayout15);
            app.UIAxes_slope2.Toolbar.Visible = 'off';
            app.UIAxes_slope2.Layout.Row = 1;
            app.UIAxes_slope2.Layout.Column = [8 9];
            app.UIAxes_slope2.Visible = 'off';

            % Create Mapdimensionrelations2Tab_2
            app.Mapdimensionrelations2Tab_2 = uitab(app.TabGroup2);
            app.Mapdimensionrelations2Tab_2.Title = 'Map-dimension relations 2';

            % Create map_dim_relations2
            app.map_dim_relations2 = uipanel(app.Mapdimensionrelations2Tab_2);
            app.map_dim_relations2.BorderType = 'none';
            app.map_dim_relations2.Position = [1 1 1340 201];

            % Create AfferentconvergenceTab_2
            app.AfferentconvergenceTab_2 = uitab(app.TabGroup2);
            app.AfferentconvergenceTab_2.Title = 'Afferent convergence';

            % Create Panel_aff_convergence
            app.Panel_aff_convergence = uipanel(app.AfferentconvergenceTab_2);
            app.Panel_aff_convergence.BorderType = 'none';
            app.Panel_aff_convergence.Position = [1 1 1340 201];

            % Create AfferentaxonarborsTab
            app.AfferentaxonarborsTab = uitab(app.TabGroup2);
            app.AfferentaxonarborsTab.Title = 'Afferent axon arbors';

            % Create Panel_aff_axon_arbors
            app.Panel_aff_axon_arbors = uipanel(app.AfferentaxonarborsTab);
            app.Panel_aff_axon_arbors.BorderType = 'none';
            app.Panel_aff_axon_arbors.Position = [1 1 1340 201];

            % Create VisualdeprivationTab
            app.VisualdeprivationTab = uitab(app.TabGroup2);
            app.VisualdeprivationTab.Title = 'Visual deprivation';

            % Create GridLayout16
            app.GridLayout16 = uigridlayout(app.VisualdeprivationTab);
            app.GridLayout16.ColumnWidth = {'1x', '1x', '1x', '1.5x', '1x', '1x'};
            app.GridLayout16.RowHeight = {'1x'};

            % Create UIAxes_orimap_contour
            app.UIAxes_orimap_contour = uiaxes(app.GridLayout16);
            app.UIAxes_orimap_contour.Toolbar.Visible = 'off';
            app.UIAxes_orimap_contour.Layout.Row = 1;
            app.UIAxes_orimap_contour.Layout.Column = 1;
            app.UIAxes_orimap_contour.Visible = 'off';

            % Create UIAxes_ori_od_hist
            app.UIAxes_ori_od_hist = uiaxes(app.GridLayout16);
            app.UIAxes_ori_od_hist.Toolbar.Visible = 'off';
            app.UIAxes_ori_od_hist.Layout.Row = 1;
            app.UIAxes_ori_od_hist.Layout.Column = 2;
            app.UIAxes_ori_od_hist.Visible = 'off';

            % Create UIAxes_orimap_biased
            app.UIAxes_orimap_biased = uiaxes(app.GridLayout16);
            app.UIAxes_orimap_biased.Layout.Row = 1;
            app.UIAxes_orimap_biased.Layout.Column = 5;
            app.UIAxes_orimap_biased.Visible = 'off';

            % Create UIAxes_od_biased
            app.UIAxes_od_biased = uiaxes(app.GridLayout16);
            app.UIAxes_od_biased.Layout.Row = 1;
            app.UIAxes_od_biased.Layout.Column = 6;
            app.UIAxes_od_biased.Visible = 'off';

            % Create UIAxes_ori_od_polar_hist
            app.UIAxes_ori_od_polar_hist = uiaxes(app.GridLayout16);
            app.UIAxes_ori_od_polar_hist.Toolbar.Visible = 'off';
            app.UIAxes_ori_od_polar_hist.Layout.Row = 1;
            app.UIAxes_ori_od_polar_hist.Layout.Column = 3;
            app.UIAxes_ori_od_polar_hist.Visible = 'off';

            % Create GridLayout_ori_biased
            app.GridLayout_ori_biased = uigridlayout(app.GridLayout16);
            app.GridLayout_ori_biased.ColumnWidth = {'3.5x', '1x'};
            app.GridLayout_ori_biased.RowHeight = {'1x', '1x', '1x', '1x'};
            app.GridLayout_ori_biased.Layout.Row = 1;
            app.GridLayout_ori_biased.Layout.Column = 4;

            % Create OrientationbiasedEditFieldLabel
            app.OrientationbiasedEditFieldLabel = uilabel(app.GridLayout_ori_biased);
            app.OrientationbiasedEditFieldLabel.HorizontalAlignment = 'right';
            app.OrientationbiasedEditFieldLabel.Layout.Row = 1;
            app.OrientationbiasedEditFieldLabel.Layout.Column = 1;
            app.OrientationbiasedEditFieldLabel.Text = 'Orientation biased';

            % Create OrientationbiasedEditField
            app.OrientationbiasedEditField = uieditfield(app.GridLayout_ori_biased, 'numeric');
            app.OrientationbiasedEditField.Layout.Row = 1;
            app.OrientationbiasedEditField.Layout.Column = 2;
            app.OrientationbiasedEditField.Value = 45;

            % Create GenerateButton
            app.GenerateButton = uibutton(app.GridLayout_ori_biased, 'push');
            app.GenerateButton.ButtonPushedFcn = createCallbackFcn(app, @GenerateButtonPushed, true);
            app.GenerateButton.Layout.Row = 4;
            app.GenerateButton.Layout.Column = [1 2];
            app.GenerateButton.Text = 'Generate';

            % Create BiaspercentDropDownLabel
            app.BiaspercentDropDownLabel = uilabel(app.GridLayout_ori_biased);
            app.BiaspercentDropDownLabel.HorizontalAlignment = 'right';
            app.BiaspercentDropDownLabel.Layout.Row = 2;
            app.BiaspercentDropDownLabel.Layout.Column = 1;
            app.BiaspercentDropDownLabel.Text = 'Bias (percent)';

            % Create BiaspercentDropDown
            app.BiaspercentDropDown = uidropdown(app.GridLayout_ori_biased);
            app.BiaspercentDropDown.Items = {'20', '30', '40', '50'};
            app.BiaspercentDropDown.Layout.Row = 2;
            app.BiaspercentDropDown.Layout.Column = 2;
            app.BiaspercentDropDown.Value = '50';

            % Create Panel
            app.Panel = uipanel(app.GridLayout);
            app.Panel.Layout.Row = 1;
            app.Panel.Layout.Column = [4 12];

            % Create GridLayout18
            app.GridLayout18 = uigridlayout(app.Panel);
            app.GridLayout18.ColumnWidth = {'1x', '1x', '1x', '1x', '1x'};
            app.GridLayout18.RowHeight = {'1x'};

            % Create UIAxes_RGC_contra
            app.UIAxes_RGC_contra = uiaxes(app.GridLayout18);
            title(app.UIAxes_RGC_contra, 'RGC Contra')
            app.UIAxes_RGC_contra.Toolbar.Visible = 'off';
            app.UIAxes_RGC_contra.Layout.Row = 1;
            app.UIAxes_RGC_contra.Layout.Column = 1;
            app.UIAxes_RGC_contra.ButtonDownFcn = createCallbackFcn(app, @CloseButtonPushed, true);
            app.UIAxes_RGC_contra.Visible = 'off';

            % Create UIAxes_RGC_ipsi
            app.UIAxes_RGC_ipsi = uiaxes(app.GridLayout18);
            title(app.UIAxes_RGC_ipsi, 'RGC Ipsi')
            app.UIAxes_RGC_ipsi.Toolbar.Visible = 'off';
            app.UIAxes_RGC_ipsi.Layout.Row = 1;
            app.UIAxes_RGC_ipsi.Layout.Column = 2;
            app.UIAxes_RGC_ipsi.Visible = 'off';

            % Create UIAxes_aff_RF_vis_space
            app.UIAxes_aff_RF_vis_space = uiaxes(app.GridLayout18);
            zlabel(app.UIAxes_aff_RF_vis_space, 'Z')
            app.UIAxes_aff_RF_vis_space.Layout.Row = 1;
            app.UIAxes_aff_RF_vis_space.Layout.Column = 3;
            app.UIAxes_aff_RF_vis_space.Visible = 'off';

            % Create UIAxes_ASF
            app.UIAxes_ASF = uiaxes(app.GridLayout18);
            zlabel(app.UIAxes_ASF, 'Z')
            app.UIAxes_ASF.Toolbar.Visible = 'off';
            app.UIAxes_ASF.Layout.Row = 1;
            app.UIAxes_ASF.Layout.Column = 4;
            app.UIAxes_ASF.Visible = 'off';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = gui_JMA_112221a

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
