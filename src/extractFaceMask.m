function [faceMask, face, ...
          filteredEyeCandidates1, ...
          filteredEyeCandidates2, ...
          filteredFaceMaskFilled] = extractFaceMask(rgbImage, ...
                                                    foregroundMask, ...
                                                    estimatedSkinMask, ...
                                                    nonFaceMask, ...
                                                    Y)
    
    faceMask = foregroundMask & estimatedSkinMask & ~nonFaceMask;

    faceMask = imfill(faceMask, 'holes');
    faceMask = extractNLargestBlobs(faceMask, 1);
    faceMask = imdilate(faceMask, strel('disk', 8));
    faceMask = imerode(faceMask, strel('disk', 10));

    faceMask = shrinkFaceMask(faceMask);

    faceMaskRep = repmat(faceMask, [1,1,3]);
    face = rgbImage.*uint8(faceMaskRep);

    grayFace = extractGrayFace(Y, face, rgbImage, faceMaskRep);
    grayFace = im2bw(grayFace, 0.47);

    filteredFaceMask = imfilter(grayFace, fspecial('laplacian'));

    filteredEyeCandidates1 = filteredFaceMask;
    filteredFaceMask = bwareaopen(filteredFaceMask, 40);
    
    originalFilteredFaceMask = filteredFaceMask;
    originalFilteredFaceMask = imdilate(originalFilteredFaceMask, strel('disk', 1));
    originalFilteredFaceMask = imerode(originalFilteredFaceMask, strel('disk', 1));
   
    [labeledImage, numberOfBlobs] = bwlabel(filteredFaceMask);
    convexProperties = regionprops(labeledImage, 'ConvexArea');
    convexAreas = [convexProperties.ConvexArea];
    [sortedConvexAreas, sortedIndices] = sort(convexAreas, 2, 'descend');
    
    filteredFaceMask = ismember(labeledImage, sortedIndices(1));
    filteredFaceMaskFilled = imfill(filteredFaceMask, 'holes');
    filteredFaceMaskFilledEroded = imerode(filteredFaceMaskFilled, strel('disk', 1));
    filteredFaceMaskBorder = xor(filteredFaceMaskFilled, filteredFaceMaskFilledEroded);

    filteredEyeCandidates1 = filteredEyeCandidates1 - filteredFaceMask;
    filteredEyeCandidates1 = imdilate(filteredEyeCandidates1, strel('disk', 5));
    filteredEyeCandidates1 = imfill(filteredEyeCandidates1, 'holes');

    filteredEyeCandidates2 = originalFilteredFaceMask - and(originalFilteredFaceMask, filteredFaceMaskBorder); 
    filteredEyeCandidates2 = imfill(filteredEyeCandidates2, 'holes');
    
    % Removed as this part has a marginal difference and therefore
    % only constraint the algorithm more than needed.
    
    % filteredFaceMask = imdilate(filteredFaceMask, strel('disk', 30));
    % filteredFaceMask = ~xor(~filteredFaceMask, imfill(filteredFaceMask, 'holes'));
    % filteredFaceMask = extractNLargestBlobs(filteredFaceMask, 1);
    % filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
    % filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
    % filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));

    % faceMask = and(filteredFaceMask, faceMask);
    % faceMaskRep = repmat(faceMask, [1,1,3]);
    % face = rgbImage.*uint8(faceMaskRep);

end


