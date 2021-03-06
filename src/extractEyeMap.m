function finalEyeMap = extractEyeMap(faceMask, filteredEyeCandidates1, ...
                                     filteredEyeCandidates2, overSaturatedMask, ...
                                     estimatedSkinMask, filteredFaceMaskEyes, ...
                                     averageFaceColorMask, grayImage, Y, Cb, Cr, mouthMap)
    
    Y = uint16(Y);
    Cr = uint16(Cr);
    Cb = uint16(Cb);
    
    eyeMapChroma = extractEyeMapChroma(Cb, Cr);

    eyeMapLuma = extractEyeMapLuma(grayImage, Y);

    eyeMap = eyeMapChroma .* eyeMapLuma;
    eyeMap = normalizeEyeMap(eyeMap);
             
    originalEyeMap = eyeMap;
    
    eyeMap = eyeMap.*filteredEyeCandidates1;
    eyeMap = eyeMap.*filteredEyeCandidates2;
    eyeMap = eyeMap.*im2double(faceMask);
    eyeMap = eyeMap.*~overSaturatedMask;
    
    eyeMap = imfill(eyeMap, 'holes');
    eyeMap = imdilate(eyeMap, strel('disk', 4));

    estimatedSkinMask2 = imerode(estimatedSkinMask, strel('disk', 5));
    estimatedSkinMask2 = imerode(estimatedSkinMask2, strel('disk', 5));
    eyeMap = eyeMap .* ~estimatedSkinMask2;
 
    eyeMap = imfill(eyeMap, 'holes');
    eyeMap = eyeMap > 0; 
    eyeMap = imdilate(eyeMap, strel('disk', 10));

    finalEyeMap = originalEyeMap;

    finalEyeMap = finalEyeMap.*eyeMap;
    finalEyeMap = finalEyeMap.*filteredFaceMaskEyes;
    finalEyeMap = finalEyeMap.*faceMask;
    finalEyeMap = finalEyeMap.*averageFaceColorMask;
    
    tentativeMouthMap = mouthMap > 0.5 * max(mouthMap(faceMask));
    tentativeMouthMap = imdilate(tentativeMouthMap, strel('disk', 7));
    tentativeMouthMap = imfill(tentativeMouthMap, 'holes');
    tentativeMouthMap = extractNLargestBlobs(tentativeMouthMap, 1);
    
    finalEyeMap = finalEyeMap .* ~tentativeMouthMap;

end


function eyeMapChroma = extractEyeMapChroma(Cb, Cr)
    CbSquared = Cb.^2;
        
    CrInvertedSquared = (1-Cr).^2;
    CbDividedByCr = Cb./Cr;
    
    eyeMapChroma = (1/3)*(CbSquared + ...
                          CrInvertedSquared + ...
                          CbDividedByCr);      
                      
    eyeMapChroma = histeq(eyeMapChroma);
    eyeMapChroma = im2double(eyeMapChroma);
end


function eyeMapLuma = extractEyeMapLuma(grayImage, Y)
    eyeMapLuma = uint16(zeros(size(grayImage)));
    for n=1 : 10
        kernel = strel('disk', n);
        dilation = imdilate(Y, kernel);
        erosion = imerode(Y, kernel);
        eyeMapLuma = eyeMapLuma + (dilation./(n + erosion));
    end
    eyeMapLuma = im2double(eyeMapLuma);
end


function eyeMap = normalizeEyeMap(eyeMap) 
    eyeMap = im2double(eyeMap);
    eyeMap = 255 * ((eyeMap - min(eyeMap(:))) ...
                 ./ (max(eyeMap(:)) - min(eyeMap(:))));
end