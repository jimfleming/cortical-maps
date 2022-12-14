function updatRetinaUI(app, appdata, eye_input)
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
    xlim(ax, [min(col_aff_all(:))-10 max(col_aff_all(:))+10])
    ylim(ax, [min(row_aff_all(:))-10 max(row_aff_all(:))+10])

    max_number_overlap_rf = max(sum_aff(:)); 
    title(ax, sprintf('Afferent sampling \n density (%.0f RF per point)', max_number_overlap_rf))

    setAxes(app, ax);
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    pause(0.001)
end

function [appdata] = show_maps_segregation(app, appdata)
    ONOFF_smoothed   = appdata.ONOFF_smoothed;
    ODCrtxPlt_smooth = appdata.ODCrtxPlt_smooth;
    map_retinotopy   = appdata.RetinotopyONOFFSortedPlot;
    cmap_map_retinotopy = app.appdata.cmap;

    onoff_od_colormap = [255 100 100; 255 0 0; 0 175 255; 0 0 255]./255;
    appdata.map1 = onoff_od_colormap;

    width_cortex = size(ODCrtxPlt_smooth, 1);
    %% figure
    axes_od = app.UIAxes_Segregation_OD ;
    axes_od_onoff = app.UIAxes_OD_ONOFF ;
    axes_retinotopy = app.UIAxes_Aff_retinotopy ;

    cla(axes_od), imagesc(ODCrtxPlt_smooth>.5, 'Parent', axes_od)
    hold(axes_od, 'on'),
    contour(axes_od, ODCrtxPlt_smooth, app.appdata.od_contour_levels, 'k', 'LineWidth', 1);
    caxis(axes_od, [0 1])
    colormap(axes_od, 'gray'), axis(axes_od, 'square')
    app.title_segregation_od.Visible = 'on';
    xlim(axes_od, [1 width_cortex]),
    ylim(axes_od, [1 width_cortex])
    setAxes(app, axes_od)

    axes_od_onoff.Visible = 'on';
    z3 = (double((ODCrtxPlt_smooth>0.5))+1)*1 ;
    z4 = (double(ONOFF_smoothed<.5)+0)*2 ;
    imagesc(z3+z4, 'Parent', axes_od_onoff)
    colormap(axes_od_onoff, onoff_od_colormap);  axis(axes_od_onoff, 'square'),
    hold(axes_od_onoff, 'on')
    app.title_eye_polarity.Visible = 'on';
    hold(axes_od_onoff, 'on'),
    contour(axes_od_onoff, ODCrtxPlt_smooth, app.appdata.od_contour_levels, 'k', 'LineWidth', 1);
    xlim(axes_od_onoff, [1 width_cortex]), ylim(axes_od, [1 width_cortex])
    setAxes(app, axes_od_onoff)

    axes_retinotopy.Visible = 'on';
    cla(axes_retinotopy)
    imshow(map_retinotopy, cmap_map_retinotopy, 'Parent', axes_retinotopy)
    axis(axes_retinotopy, 'square')
    app.title_retinotopy.Visible = 'on';
    xlim(axes_retinotopy, [1 width_cortex]), ylim(axes_retinotopy, [1 width_cortex])
    setAxes(app, axes_retinotopy)
end

function update_primord_orimap_UI(app, appdata)
    app.title_primord_orimap.Visible = 'on';

    primord_orimap_contra = appdata.data_primord_contra.OriPreferred;
    primord_orimap_ipsi = appdata.data_primord_ipsi.OriPreferred;

    axes1 = app.UIAxes_primord_orimap_contra;
    cla(axes1)
    imshow(primord_orimap_contra, hsv(180), 'Parent', axes1); hold(axes1, 'on');
    contour( axes1, appdata.ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
    xlabel(axes1, 'contra') % , 'fontsize', 14)
    setAxes(app, axes1);
    axes1.XColor = appdata.contra_font_color;
    axes1.YColor = appdata.contra_font_color;

    axes2 = app.UIAxes_primord_orimap_ipsi;
    cla(axes2)
    imshow(primord_orimap_ipsi, hsv(180), 'Parent', axes2); hold(axes2, 'on');
    contour(axes2, appdata.ODCrtxPlt_smooth, 1, 'k', 'LineWidth', 1);
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
    xlabel(ax1, 'contra')
    contour(ax1, appdata.ODCrtxPlt_interpolated, app.appdata.od_contour, 'k', 'LineWidth', 1);% app.od_contour_w);
    setAxes(app, ax1);

    ax2 = app.UIAxes_mature_orimap_ipsi;
    imshow(mature_orimap_ipsi, hsv(180), 'Parent', ax2); hold(ax2,'on');
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
    
    app.Panel_aff_convergence.Visible = 'on'; 
    t2 = tiledlayout(app.Panel_aff_convergence, 2, 11);

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
    angles     = OriBinRFReference{cortex_row, cortex_col};
    ori_tuning = OriHistRFReference{cortex_row, cortex_col};
    ori_tuning_odi = ori_tuning * ODI(cortex_row, cortex_col);
    polarS3(ax, angles, ori_tuning_odi', axis_color);
    pause(0.001)
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

    hold(ax,'on')
    rectangle(ax, 'Position', [hh-CenterRadiusX/2, hh-CenterRadiusY/2, CenterRadiusX, CenterRadiusY],'EdgeColor', 'k', 'LineWidth', 2, 'curvature', [1 1])
    rectangle(ax, 'Position', [hh-SurroundRadiusX/2, hh-SurroundRadiusY/2, SurroundRadiusX, SurroundRadiusY],'EdgeColor', 'k', 'LineWidth', 2, 'curvature', [1 1])
    set(ax, 'ydir', 'reverse')

    title(ax, sprintf('Afferent sorting \n filter'))
    axis(ax, 'square')
    xlim(ax, [0 CrtxLength])
    ylim(ax, [0 CrtxLength])

    setAxes(app, ax);
    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];
    pause(0.001)
end

function StartButtonPushed(app, event)
    % Loading user's main parameters
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
    ADF = str2double(app.AffdenfactorDropDown.Value);
    app.appdata.retina.aff_density_factor = ADF; 
    app.appdata.pix2mic = 22;

    % Loading user's cortical parameters
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
    app.appdata.pix2deg = 40;

    % Loading user's afferent sorting parameters
    app.appdata.NSortOD = 10;    % app.ODsortingiterationsEditField.Value;               % Number of iteration for OD sorting
    app.appdata.NSortONOFF = 10; % app.ONOFFsortingiterationsEditField.Value;            % Number of iteration for ON/OFF sorting

    % Sorting filter parameters come from the RGC distribution
    % CenterRadiusX is measured based on RGC distance in make_retina_app function 
    aff_sampling_density        = app.appdata.CenterRadiusX;
    app.appdata.aff_sampling_density = aff_sampling_density; 
    
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

    app.appdata.r_covered_aff = 20; %20; % average radius covered by RFs in cortex in pixel

    rf_sd_RGC = ADF * 5; % default val : 5;
    app.appdata.RFONCenter = rf_sd_RGC; % 5;                        % Receptive Field Center
    app.appdata.RFONSurround = round(1.1*app.appdata.RFONCenter);  % Receptive Field Surround
    app.appdata.RFOFFCenter = app.appdata.RFONCenter;
    app.appdata.RFOFFSurround = app.appdata.RFONSurround;

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
    [app.appdata] = show_maps_segregation(app, app.appdata);
    update_asf_figure(app, app.appdata)

    % labeling islands
    [app.appdata.ONOFFODLabelSorted, app.appdata.CenterPinwheelSorted, app.appdata.NumPinwheel, app.appdata.center_island_image] = label_ononff_od_island_gui(app.appdata, debug, 0);

    % synaptic weight based on islands
    app.appdata.weight_map_island = ones(size(app.appdata.ONOFFODLabelSorted)); % it should be tested whether adding weights here helps bring the pw closer to the center of islands

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

        app.appdata.eyepref = 'ipsi';   % default ipsi
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
