function [mouthX, mouthY] = extractMouth(faceMask, face, nonMouthMask)

    ycbcrFace = rgb2ycbcr(face);
    Cb = double(ycbcrFace(:,:,2));
    Cr = double(ycbcrFace(:,:,3));
    
    n = length(nonzeros(faceMask));
    
    CrSquared = Cr.^2;
    numerator = (1 / n) * sum(CrSquared(:));
    
    CrDividedByCb = Cr./Cb;
    denumerator = (1 / n) * sum(CrDividedByCb(:));
    
    nn = 0.95 * (numerator / denumerator);
    
    mouthMap = CrSquared.*(CrSquared - nn*CrDividedByCb).^2;
    
    mouthMap = im2double(mouthMap);
    mouthMap = 255 * ((mouthMap - min(mouthMap(:))) ...
                 ./ (max(mouthMap(:)) - min(mouthMap(:))));
    mouthMap = uint8(mouthMap);
    
    mouthMap = mouthMap .* uint8(faceMask);
%     figure; imshow(mouthMap); title('mouthMap'); pause;

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