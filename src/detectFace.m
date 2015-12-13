function output = detectFace(rgbImage)

    rgbImage = padarray(rgbImage, [2 2], 'both');
    
    grayImage = rgb2gray(rgbImage);
    
    hsvImage = rgb2hsv(rgbImage);
    
    ycbcrImage = rgb2ycbcr(rgbImage);

    Y = ycbcrImage(:,:,1);
    Cb = ycbcrImage(:,:,2);
    Cr = ycbcrImage(:,:,3);

    [estimatedSkinMask] = extractSkinColorFromCbCr(Cb, Cr);

    overSaturatedMask = extractOverSaturated(grayImage);

    foregroundMask = extractForeground(grayImage, overSaturatedMask);
    
    nonFaceMask = extractNonFaceMask(Cb, Cr, hsvImage);
    
    try
        [faceMask, ...
         face, ... 
         filteredEyeCandidates1, ...
         filteredEyeCandidates2, ...
         filteredFaceMaskEyes] = extractFaceMask(rgbImage, ...
                                                 foregroundMask, ....
                                                 estimatedSkinMask, ...
                                                 nonFaceMask, ...
                                                 Y);                                 
    catch exception
        faceMask = foregroundMask & estimatedSkinMask & ~nonFaceMask;
        faceMask = imfill(faceMask, 'holes');
        faceMask = extractNLargestBlobs(faceMask, 1);
        faceMask = imdilate(faceMask, strel('disk', 8));
        faceMask = imerode(faceMask, strel('disk', 10));
        faceMask = shrinkFaceMask(faceMask);
        faceMaskRep = repmat(faceMask, [1,1,3]);
        face = rgbImage.*uint8(faceMaskRep);
        [row, col] = find(face(:,:,1) ~= 0);
        output = rgbImage(min(row):max(row), min(col):max(col), :);
        return
    end
    
    try
        averageFaceColorMask = extractAverageFaceColor(rgbImage, faceMask);

        mouthMap = extractMouthMap(faceMask, face);

        eyeMap = extractEyeMap(faceMask, filteredEyeCandidates1, ...
                               filteredEyeCandidates2, overSaturatedMask, ...
                               estimatedSkinMask, filteredFaceMaskEyes, ...
                               averageFaceColorMask, grayImage, ...
                               Y, Cb, Cr, mouthMap);

        [leftEyeCenterX, leftEyeCenterY, leftEyeRadius, ...
         rightEyeCenterX, rightEyeCenterY, rightEyeRadius] = extractEyes(eyeMap);

        [leftEyeCenterX, leftEyeCenterY, ...
         rightEyeCenterX, rightEyeCenterY, ...
         rgbImage, face, mouthMap] = rotateFace(leftEyeCenterX, leftEyeCenterY, ...
                                                rightEyeCenterX, rightEyeCenterY, ...
                                                rgbImage, face, mouthMap);

        eyesMeanY = (leftEyeCenterY + rightEyeCenterY) / 2;
        eyeDistance = rightEyeCenterX - leftEyeCenterX;
        offsetX = eyeDistance * 0.3;

        nonMouthMask = extractNonMouthMask(face, eyeDistance, eyesMeanY, ...
                                           leftEyeCenterX, rightEyeCenterX, ...
                                           offsetX);

        [mouthX, mouthY] = extractMouth(mouthMap, nonMouthMask);

        eyeMouthDistance = mouthY - eyesMeanY;
        offsetY = eyeMouthDistance * 0.2; 

        output = rgbImage(round(eyesMeanY-offsetY) : round(mouthY+offsetY), ...
                          round(leftEyeCenterX-offsetX) : round(rightEyeCenterX+offsetX), :);

    catch exception
        [row, col] = find(face(:,:,1) ~= 0);
        output = rgbImage(min(row):max(row), min(col):max(col), :);
    end
    
    figure; imshow(output); title('output'); 
    
end