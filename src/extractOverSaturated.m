function overSaturatedMask = extractOverSaturated(grayImage)
    
    overSaturated = imfilter(grayImage, fspecial('laplacian'));
    overSaturatedMask = ~imdilate(overSaturated, strel('disk', 2));

end
