function output = detectFace(rgbImage)

    rgbImage = padarray(rgbImage, [2 2], 'both');
%     figure; imshow(rgbImage); title('rgbImage');    

    grayImage = rgb2gray(rgbImage);
%     figure; imshow(grayImage); title('grayImage'); pause;

    hsvImage = rgb2hsv(rgbImage);
%     figure; imshow(hsvImage); title('hsvImage'); pause;

    [estimatedSkinMask, Y, Cb, Cr] = extractSkinColorFromYCbCr(rgbImage);
%     figure; imshow(estimatedSkinMask); title('estimatedSkinMask'); pause;   
    
    overSaturatedMask = extractOverSaturated(grayImage);
%     figure; imshow(overSaturatedMask); title('overSaturatedMask'); pause;   
    
    foregroundMask = extractForeground(grayImage, overSaturatedMask, rgbImage);
%     figure; imshow(foregroundMask); title('foregroundMask'); pause; 
    
    nonFaceMask = extractNonFaceMask(Cb, Cr, hsvImage);
%     figure; imshow(nonFaceMask); title('nonFaceMask'); pause;

    [faceMask, face, filteredFaceMaskCopy, filteredFaceMaskCopy2, ...
     filteredFaceMaskEyes] = extractFaceMask(rgbImage, foregroundMask, ....
                                             estimatedSkinMask, nonFaceMask);
    
    averageFaceColorMask = extractAverageFaceColor(rgbImage, faceMask);
%     figure; imshow(averageFaceColorMask); title('averageFaceColorMask'); pause;

    finalEyeMap = extractEyeMap(faceMask, filteredFaceMaskCopy, ...
                                filteredFaceMaskCopy2, overSaturatedMask, ...
                                estimatedSkinMask, filteredFaceMaskEyes, ...
                                averageFaceColorMask, grayImage, Y, Cb, Cr);
  %     figure; imshow(finalEyeMap); title('finalEyeMap'); pause;
    
    [leftEyeCenterX, leftEyeCenterY, leftEyeRadius, ...
     rightEyeCenterX, rightEyeCenterY, rightEyeRadius] = extractEyes(finalEyeMap);
    
    [leftEyeCenterX, leftEyeCenterY, ...
     rightEyeCenterX, rightEyeCenterY, ...
     rgbImage, face, faceMask] = rotateFace(leftEyeCenterX, leftEyeCenterY, ...
                                            rightEyeCenterX, rightEyeCenterY, ...
                                            rgbImage, face, faceMask);
    
    eyesMeanY = (leftEyeCenterY + rightEyeCenterY) / 2;
    eyeDistance = rightEyeCenterX - leftEyeCenterX;
    offsetX = eyeDistance * 0.3; % 0.3
    
    nonMouthMask = extractNonMouthMask(face, eyeDistance, eyesMeanY, ...
                                       leftEyeCenterX, rightEyeCenterX, ...
                                       offsetX);
    
    [mouthX, mouthY] = extractMouth(faceMask, face, nonMouthMask);

    eyeMouthDistance = mouthY - eyesMeanY;
    offsetY = eyeMouthDistance * 0.2; % 0.2
    
    output = rgbImage(round(eyesMeanY-offsetY) : round(mouthY+offsetY), ...
                      round(leftEyeCenterX-offsetX) : round(rightEyeCenterX+offsetX), :);
    
%     figure; imshow(output); title('output'); %pause;
    
end