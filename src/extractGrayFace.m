function grayFace = extractGrayFace(Y, face, rgbImage, faceMaskRep)

    rgbImage2 = rgbImage;
    rgbImage(:,:,1) = rgbImage(:,:,1) .* ((255 - Y) ./ 255);
    rgbImage(:,:,2) = rgbImage(:,:,2) .* ((255 - Y) ./ 255);
    rgbImage(:,:,3) = rgbImage(:,:,3) .* ((255 - Y) ./ 255);

    dark = rgbImage .* uint8(faceMaskRep);
    bright = (rgbImage2 - dark) .*  uint8(faceMaskRep);

    bRmask = find(bright(:,:,1)>0);
    bGmask = find(bright(:,:,2)>0);
    bBmask = find(bright(:,:,3)>0);
    
    bRm = mean( bright(bRmask) );
    bGm = mean( bright(bGmask) );
    bBm = mean( bright(bBmask) );
    
    dRmask = find(dark(:,:,1)>0);
    dGmask = find(dark(:,:,2)>0);
    dBmask = find(dark(:,:,3)>0);
    
    dRm = mean( dark(dRmask) );
    dGm = mean( dark(dGmask) );
    dBm = mean( dark(dBmask) );

    redDiff = bRm - dRm;
    greenDiff = bGm - dGm;
    blueDiff = bBm - dBm;
    
    sumRscale= (sum(bRmask(:)) / sum(dRmask(:)));
    sumGscale = (sum(bGmask(:)) / sum(dGmask(:)));
    sumBscale = (sum(bBmask(:)) / sum(dBmask(:)));
    
    twoFace = (sumRscale < 0.7) && ...
              (sumGscale < 0.7) && ...
              (sumBscale < 0.7) && ...
              redDiff > 75 && redDiff < 120 && ...
              greenDiff > 75 && greenDiff < 120 && ...
              blueDiff > 75 && blueDiff < 120;
        
%     figure; imshow(bright); title('bright'); pause;
%     figure; imshow(dark); title('dark'); pause;
       
    redCont = uint8(zeros(size(dark(:,:,1))));
    redCont(dRmask) = 0.8 * redDiff;
    
    greenCont = uint8(zeros(size(dark(:,:,2))));
    greenCont(dGmask) = 0.8 * greenDiff;
    
    blueCont = uint8(zeros(size(dark(:,:,3))));
    blueCont(dBmask) = 0.8 *  blueDiff;
      
    dark2 = cat( 3, ( bright(:,:,1) + dark(:,:,1) + redCont ) , ...
                      bright(:,:,2) + dark(:,:,2) + greenCont, ...
                      bright(:,:,2) + dark(:,:,3) + blueCont);
            
    if twoFace
        grayFace = rgb2gray(dark2);
    else
        grayFace = imadjust(rgb2gray(face));
    end
    
end