function [mouthX, mouthY] = extractMouth(mouthMap, nonMouthMask)

%     figure; imshow(mouthMap); title('mouthMap'); pause;
%     figure; imshow(~nonMouthMask); pause; %title('~nonMouthMask'); pause;
    
    mouthMap = mouthMap > 0.5 * max(mouthMap(~nonMouthMask));

    mouthMap = imdilate(mouthMap, strel('disk', 4));
    mouthMap = imfill(mouthMap, 'holes');
    mouthMap = imerode(mouthMap, strel('disk', 4));
    
    mouthMask = and(mouthMap, ~nonMouthMask);
    
    try
    mouthMask = ExtractNLargestBlobs(mouthMask, 1);
    catch
        warning('extractMouthMap.m failed. At ExtractNLargestBlobs.');
        mouthMask = -1;
        return;
    end
    
    mouthBlobs = regionprops(mouthMask, mouthMap, {'Centroid'});

    mouthX = mouthBlobs(1).Centroid(1);
    mouthY = mouthBlobs(1).Centroid(2);
    
end