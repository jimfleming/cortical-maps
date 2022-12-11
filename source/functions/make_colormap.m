function [cmap,RetinotopyRFspace_plot] = make_colormap(CrtxGridXX,debug)

    nRGB1 = size(CrtxGridXX,1); 
    NColor = nRGB1; 
    RetinotopyRFspace_plot = zeros(nRGB1,nRGB1);
    RetinotopyRFspace_plot(:) = 1 : nRGB1 * nRGB1 ; 
    [RedGrid,~] = meshgrid( linspace(0,255,NColor) , linspace(0,0,NColor) ); 
    [~,GreenGrid] = meshgrid( linspace(230,230,NColor) , linspace(230,51,NColor) ); 
    [BlueGrid,~] = meshgrid( linspace(0,255,NColor) , linspace(0,0,NColor) ); 
    cmap = [RedGrid(:) GreenGrid(:) BlueGrid(:)]/255;
    
    if debug == 1 
       figure, imshow(RetinotopyRFspace_plot,cmap)
    end
end 