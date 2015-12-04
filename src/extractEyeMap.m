function finalEyeMap = extractEyeMap(faceMask, filteredFaceMaskCopy, ...
                                     filteredFaceMaskCopy2, overSaturatedMask, ...
                                     estimatedSkinMask, filteredFaceMaskEyes, ...
                                     averageFaceColorMask, grayImage, Y, Cb, Cr)
    
    Y = uint16(Y);
    Cr = uint16(Cr); % Cb also works
    Cb = uint16(Cb);
    
    eyeMapChroma = extractEyeMapChroma(Cb, Cr);
    
    eyeMapLuma = extractEyeMapLuma(grayImage, Y);
    
    eyeMap = eyeMapChroma .* eyeMapLuma;
    
    eyeMap = normalizeEyeMap(eyeMap);
             
    originalEyeMap = eyeMap;
%     figure; imshow(eyeMap, []); title('originalEyeMap'); pause;

    eyeMap = eyeMap.*filteredFaceMaskCopy;
%     figure; imshow(filteredFaceMaskCopy); title('filteredFaceMaskCopy'); pause;
%     figure; imshow(eyeMap); title('eyeMap.*filteredFaceMaskCopy'); pause;
    
    eyeMap = eyeMap.*filteredFaceMaskCopy2;
%     figure; imshow(filteredFaceMaskCopy2); title('filteredFaceMaskCopy2'); pause;
%     figure; imshow(eyeMap); title('eyeMap.*filteredFaceMaskCopy2'); pause;
    
    eyeMap = eyeMap.*im2double(faceMask);
%     figure; imshow(eyeMap); title('eyeMap.*im2double(faceMask)'); pause;
    


    eyeMap = eyeMap.*~overSaturatedMask;
%     figure; imshow(eyeMap); title('eyeMap'); pause;
    eyeMap = imfill(eyeMap, 'holes');
    eyeMap = imdilate(eyeMap, strel('disk', 4));

%     figure; imshow(eyeMap); title('before final Eyemap'); pause;
    
    

    estimatedSkinMask2 = imerode(estimatedSkinMask, strel('disk', 5));
    estimatedSkinMask2 = imerode(estimatedSkinMask2, strel('disk', 5));
    eyeMap = eyeMap .* ~estimatedSkinMask2;
    
%     figure; imshow(eyeMap); title('eyeMap'); pause;

    eyeMap = imfill(eyeMap, 'holes');
    eyeMap = eyeMap > 0; 
    eyeMap = imdilate(eyeMap, strel('disk', 10));

    finalEyeMap = originalEyeMap;
    finalEyeMap = finalEyeMap.*eyeMap;
    finalEyeMap = finalEyeMap.*filteredFaceMaskEyes;
    finalEyeMap = finalEyeMap.*faceMask;
    finalEyeMap = finalEyeMap.*averageFaceColorMask;

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
%         figure; imshow(eyeMapLuma); title('eyeMapLuma'); pause;
    end

    eyeMapLuma = im2double(eyeMapLuma);
end


function eyeMap = normalizeEyeMap(eyeMap) 
    eyeMap = im2double(eyeMap);
    eyeMap = 255 * ((eyeMap - min(eyeMap(:))) ...
                 ./ (max(eyeMap(:)) - min(eyeMap(:))));
end