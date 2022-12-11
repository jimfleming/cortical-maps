function [ONOFFODLabelSorted,CenterPinwheelSorted,NumPinwheel,CenterPWImage] = label_ononff_od_island2(ODCrtxPlt,ONOFFCrtxPlt,debug, show_fig)

    % Labeling each ONOFF OD feature separately
    FeatureRemoveCriteria = 8 ; % in pixel (islands with area smaller than this number are removed)
    disk_openning_size = 0; 

    [ONOFFODLabel] = label_onoff_unsorted(ODCrtxPlt,ONOFFCrtxPlt,disk_openning_size);
    [FeatureCenter,FeatureArea,DistToCenter,CenterPWImage] = calculate_center_island(ONOFFODLabel,FeatureRemoveCriteria,debug);

    % Resort the pinwheels index based on the distance to the center of image
    FeatureIncluded = FeatureArea > FeatureRemoveCriteria;% in pixel, removing small islands 
    NumPinwheel = sum(FeatureIncluded(:));
    [~,indPinwheelSorted] = sort(DistToCenter);
    ONOFFODLabelSorted = zeros(size(ONOFFODLabel,1) , size(ONOFFODLabel,2));
    CenterPinwheelSorted{NumPinwheel} = [];
    cc = 0; 
    for bb = 1 : length(indPinwheelSorted) %  NumPinwheel
        iPw = indPinwheelSorted(bb); 
        if FeatureArea(iPw) > FeatureRemoveCriteria
            cc = cc + 1 ; 
            ONOFFODLabelSorted = ONOFFODLabelSorted + (ONOFFODLabel == iPw) .* cc;
            CenterPinwheelSorted{cc} =  FeatureCenter{iPw};
        end
    end
    
    % Plot ON/OFF OD islands 
    if show_fig == 1
        figure,imagesc(ONOFFODLabelSorted),colorbar,colormap('jet')
        hold on
        for ii = 1 : NumPinwheel
            plot(CenterPinwheelSorted{ii}(2),CenterPinwheelSorted{ii}(1),'ko')
        end
        axis square
        title('sorted island')
    end

end 


function [ONOFFODLabel] = label_onoff_unsorted(ODCrtxPlt,ONOFFCrtxPlt,disk_openning_size)
    SEOpenning = strel('disk',disk_openning_size); 
    BWNumConnection = 4; 

    ONOFFODTemp = ONOFFCrtxPlt*2 + ODCrtxPlt ;
    IpsiOFFOpened = imopen(ONOFFODTemp == -3,SEOpenning);
    [IpsiOFFLabeled, IpsiOFFNum] = bwlabel(IpsiOFFOpened,BWNumConnection);%   figure,imagesc(IpsiOFFLabeled),colorbar
    ContraOFFOpened = imopen(ONOFFODTemp == -1,SEOpenning);
    [ContraOFFLabeled, ContraOFFNum] = bwlabel(ContraOFFOpened,BWNumConnection);%   figure,imagesc(ContraOFFLabeled),colorbar
    ContraONOpened = imopen(ONOFFODTemp == 3,SEOpenning);
    [ContraONLabeled, ContraONNum] = bwlabel(ContraONOpened,BWNumConnection); 
    IpsiONOpened = imopen(ONOFFODTemp == 1,SEOpenning);
    [IpsiONLabeled, ~] = bwlabel(IpsiONOpened,BWNumConnection);

    %   Labelling again giving each ON OFF Ipsi Contra a number
    ONOFFODLabel1 = IpsiOFFLabeled ;  %   figure,imagesc(ONOFFODLabel1)
    ONOFFODLabel2 = ContraOFFLabeled + (ContraOFFLabeled ~= 0)*(IpsiOFFNum) ; %   figure,imagesc(ONOFFODLabel2),colorbar
    ONOFFODLabel3 = ContraONLabeled + (ContraONLabeled ~= 0)*(IpsiOFFNum+ContraOFFNum) ; %   figure,imagesc(ONOFFODLabel3)
    ONOFFODLabel4 = IpsiONLabeled + (IpsiONLabeled ~= 0)*(IpsiOFFNum+ContraOFFNum+ContraONNum) ; %   figure,imagesc(ONOFFODLabel4)
    ONOFFODLabel = ONOFFODLabel1 + ONOFFODLabel2 + ONOFFODLabel3 + ONOFFODLabel4;% figure,imagesc(ONOFFODLabel),colorbar,colormap('jet')
end 
    

function [FeatureCenter,FeatureArea,DistToCenter,CenterPWImage] = calculate_center_island(ONOFFODLabel,RemoveCriteria,debug)
    % Finding the Center of each ONOFF OD feature
    % Calculating distance from center of Islands (for Spatial Frequency measurement)

    ONOFFODNum = max(ONOFFODLabel(:));
    FeatureCenter{ONOFFODNum} = [];
    FeatureArea = zeros(ONOFFODNum,1);
    DistToCenter = zeros(ONOFFODNum,1);
    RowHalfSize = size(ONOFFODLabel,1)/2;
    ColHalfSize = size(ONOFFODLabel,2)/2;
    CenterPWImage = zeros(size(ONOFFODLabel)); 
    for ii = 1 : ONOFFODNum
        Feature = ONOFFODLabel == ii; 
        Center = regionprops(Feature,'centroid');

        % The result of regionprops could be out of the region so the closest point to the central line (skeletion) is found
        ColCenterMidLine = Center.Centroid(1);
        RowCenterMidLine = Center.Centroid(2);
        FeatureSkel = bwmorph(Feature,'skel',Inf);
        [r1,c1] = find(FeatureSkel);
        dist = (c1-ColCenterMidLine).^2 + (r1-RowCenterMidLine).^2;
        [~,aa] = min(dist);
        FeatureCenter{ii} = [r1(aa) c1(aa)];
        FeatureArea(ii) = sum(Feature(:));
        DistToCenter(ii) = sqrt((c1(aa)-ColHalfSize).^2 + (r1(aa)-RowHalfSize).^2);

        if FeatureArea(ii) > RemoveCriteria
            CenterPWImage(r1(aa),c1(aa)) = 1; %   figure,imagesc(CenterPWImage)
        end 
    end

    if debug == 1 
        figure,imagesc(ONOFFODLabel),colorbar,colormap('jet')
        hold on
        for ii = 1 : ONOFFODNum
            plot(FeatureCenter{ii}(2),FeatureCenter{ii}(1),'o')
        end
        title('Unsorted islands')
    end 
end 
