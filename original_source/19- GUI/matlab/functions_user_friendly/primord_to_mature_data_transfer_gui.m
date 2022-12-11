function appdata = primord_to_mature_data_transfer_gui(appdata)

n_ori_smooth = appdata.n_ori_smooth; 
    %%  map analysis (Map interpolation)
appdata.data_contra_mature.SF_map   = appdata.data_primord_contra.SFPref;
appdata.data_contra_mature.SF50_map = appdata.data_primord_contra.SF50Pref;
appdata.data_contra_mature.CV_map   = appdata.data_primord_contra.CV;
appdata.data_contra_mature.LPI_map  = appdata.data_primord_contra.LPI;
appdata.data_contra_mature.LHI_map  = appdata.data_primord_contra.LHI_map;

appdata.data_contra_mature.max_on_response  = appdata.data_primord_contra.max_on_response;
appdata.data_contra_mature.max_off_response  = appdata.data_primord_contra.max_off_response;
appdata.data_contra_mature.max_on_response_norm  = appdata.data_primord_contra.max_on_response_norm;
appdata.data_contra_mature.max_off_response_norm  = appdata.data_primord_contra.max_off_response_norm;

appdata.data_contra_mature.sf_map_interpolated   = imresize(appdata.data_primord_contra.SFPref, n_ori_smooth);
appdata.data_contra_mature.sf50_map_interpolated = imresize(appdata.data_primord_contra.SF50Pref, n_ori_smooth);
appdata.data_contra_mature.cv_map_interpolated   = imresize(appdata.data_primord_contra.CV, n_ori_smooth);
appdata.data_contra_mature.LPI_intepolated       = imresize(appdata.data_primord_contra.LPI, n_ori_smooth);
appdata.data_contra_mature.LHI_interpolated      = imresize(appdata.data_primord_contra.LHI_map, n_ori_smooth);

%%
appdata.data_ipsi_mature.SF_map   = appdata.data_primord_ipsi.SFPref;
appdata.data_ipsi_mature.SF50_map = appdata.data_primord_ipsi.SF50Pref;
appdata.data_ipsi_mature.CV_map   = appdata.data_primord_ipsi.CV;
appdata.data_ipsi_mature.LPI_map  = appdata.data_primord_ipsi.LPI;
appdata.data_ipsi_mature.LHI_map  = appdata.data_primord_ipsi.LHI_map;

appdata.data_ipsi_mature.max_on_response  = appdata.data_primord_ipsi.max_on_response;
appdata.data_ipsi_mature.max_off_response  = appdata.data_primord_ipsi.max_off_response;
appdata.data_ipsi_mature.max_on_response_norm  = appdata.data_primord_ipsi.max_on_response_norm;
appdata.data_ipsi_mature.max_off_response_norm  = appdata.data_primord_ipsi.max_off_response_norm;

appdata.data_ipsi_mature.sf_map_interpolated   = imresize(appdata.data_primord_ipsi.SFPref, n_ori_smooth);
appdata.data_ipsi_mature.sf50_map_interpolated = imresize(appdata.data_primord_ipsi.SF50Pref, n_ori_smooth);
appdata.data_ipsi_mature.cv_map_interpolated   = imresize(appdata.data_primord_ipsi.CV, n_ori_smooth);
appdata.data_ipsi_mature.LPI_intepolated       = imresize(appdata.data_primord_ipsi.LPI, n_ori_smooth);
appdata.data_ipsi_mature.LHI_interpolated      = imresize(appdata.data_primord_ipsi.LHI_map, n_ori_smooth);


end 
