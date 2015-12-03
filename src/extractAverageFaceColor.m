function averageFaceColorMask = extractAverageFaceColor(rgbImage, faceMask) 

    averageFaceColorMask = meanColorize(rgbImage, faceMask, 0.2);
    averageFaceColorMask = imerode(averageFaceColorMask, strel('disk', 2));
    averageFaceColorMask = imfill(averageFaceColorMask, 'holes');
%     figure; imshow(averageFaceColorMask); title('averageFaceColorMask'); 

end