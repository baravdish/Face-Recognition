function [mouthX, mouthY] = extractMouth(mouthMap, nonMouthMask)

    mouthMap = mouthMap > 0.5 * max(mouthMap(~nonMouthMask));

    mouthMap = imdilate(mouthMap, strel('disk', 4));
    mouthMap = imfill(mouthMap, 'holes');
    mouthMap = imerode(mouthMap, strel('disk', 4));
    
    mouthMask = and(mouthMap, ~nonMouthMask);
    mouthMask = extractNLargestBlobs(mouthMask, 1);
    
    mouthBlobs = regionprops(mouthMask, mouthMap, {'Centroid'});

    mouthX = mouthBlobs(1).Centroid(1);
    mouthY = mouthBlobs(1).Centroid(2);
    
end