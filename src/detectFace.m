function output = detectFace(rgbImage, originalRgbImage)
    
%     figure; imshow(rgbImage); title('rgbImage')
    
    rgbImage = padarray(rgbImage, [2 2], 'both');
%     figure; imshow(rgbImage); title('padded rgbImage');

    grayImage = rgb2gray(rgbImage);
%     figure; imshow(grayImage); title('grayImage'); pause;

    hsvImage = rgb2hsv(rgbImage);
%     figure; imshow(hsvImage); title('hsvImage'); pause;

    ycbcrImage = rgb2ycbcr(rgbImage);
%     figure; imshow(ycbcrImage); title('ycbcrImage'); pause;

    Y = ycbcrImage(:,:,1);
    Cb = ycbcrImage(:,:,2);
    Cr = ycbcrImage(:,:,3);

    [estimatedSkinMask] = extractSkinColorFromCbCr(Cb, Cr);
%     figure; imshow(estimatedSkinMask); title('estimatedSkinMask'); pause;

    overSaturatedMask = extractOverSaturated(grayImage);
%     figure; imshow(overSaturatedMask); title('overSaturatedMask'); pause;

    foregroundMask = extractForeground(grayImage, overSaturatedMask, rgbImage);
%     figure; imshow(foregroundMask); title('foregroundMask'); pause;
%     figure; imshow(~foregroundMask); title('backgroundMask'); pause;
    
    nonFaceMask = extractNonFaceMask(Cb, Cr, hsvImage);
%     figure; imshow(nonFaceMask); title('nonFaceMask'); pause;

    [faceMask, face, filteredEyeCandidates1, filteredEyeCandidates2, ...
     filteredFaceMaskEyes] = extractFaceMask(rgbImage, foregroundMask, ....
                                             estimatedSkinMask, nonFaceMask);
%     figure; imshow(face); title('face'); pause;
    
    averageFaceColorMask = extractAverageFaceColor(rgbImage, faceMask);
%     figure; imshow(averageFaceColorMask); title('averageFaceColorMask'); pause;

    eyeMap = extractEyeMap(faceMask, filteredEyeCandidates1, ...
                           filteredEyeCandidates2, overSaturatedMask, ...
                           estimatedSkinMask, filteredFaceMaskEyes, ...
                           averageFaceColorMask, grayImage, Y, Cb, Cr);
%     figure; imshow(eyeMap./max(eyeMap(:))); title('eyeMap'); pause;
    
    [leftEyeCenterX, leftEyeCenterY, leftEyeRadius, ...
     rightEyeCenterX, rightEyeCenterY, rightEyeRadius] = extractEyes(eyeMap);
    
    [leftEyeCenterX, leftEyeCenterY, ...
     rightEyeCenterX, rightEyeCenterY, ...
     rgbImage, face, faceMask] = rotateFace(leftEyeCenterX, leftEyeCenterY, ...
                                            rightEyeCenterX, rightEyeCenterY, ...
                                            rgbImage, face, faceMask, ...
                                            leftEyeRadius, rightEyeRadius);

    eyesMeanY = (leftEyeCenterY + rightEyeCenterY) / 2;
    eyeDistance = rightEyeCenterX - leftEyeCenterX;
    offsetX = eyeDistance * 0.3; % 0.3

    nonMouthMask = extractNonMouthMask(face, eyeDistance, eyesMeanY, ...
                                       leftEyeCenterX, rightEyeCenterX, ...
                                       offsetX);
%     figure; imshow(nonMouthMask); title('nonMouthMask'); pause;
%     figure; imshow(~nonMouthMask); title('~nonMouthMask'); pause;
    
    [mouthX, mouthY] = extractMouth(faceMask, face, nonMouthMask);

    eyeMouthDistance = mouthY - eyesMeanY;
    offsetY = eyeMouthDistance * 0.2; % 0.2
    
%     output = rgbImage(round(eyesMeanY-offsetY) : round(mouthY-2*offsetY), ...
%                       round(leftEyeCenterX-offsetX) : round(rightEyeCenterX+offsetX), :);
        output = rgbImage(round(eyesMeanY-offsetY) : round(mouthY+offsetY), ...
                      round(leftEyeCenterX-offsetX) : round(rightEyeCenterX+offsetX), :);

%     output = lightNorm(output);
    
    figure; imshow(output); title('output'); %pause;

end

