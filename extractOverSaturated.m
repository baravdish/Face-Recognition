function overSaturatedMask = extractOverSaturated(grayImage)
    
    overSaturated = imfilter(grayImage, fspecial('laplacian'));
    overSaturatedMask = ~imdilate(overSaturated, strel('disk', 2));

%     figure; imshow(grayImage); title('rgbImage before detracting over saturation'); pause;
%     grayImage = grayImage.*uint8(~overSaturatedMask);
%     figure;  imshow(grayImage); title('rgbImage after detracting over saturation'); pause;

end
