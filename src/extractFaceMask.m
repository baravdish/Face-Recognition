function [faceMask, face, ...
          filteredEyeCandidates1, ...
          filteredEyeCandidates2, ...
          filteredFaceMaskFilled] = extractFaceMask(rgbImage, ...
                                                    foregroundMask, ...
                                                    estimatedSkinMask, ...
                                                    nonFaceMask, ...
                                                    Y)
    
    faceMask = foregroundMask & estimatedSkinMask & ~nonFaceMask;
%     figure; imshow(faceMask); title('foregroundMask & estimatedSkinMask & ~nonFaceMask'); pause;
   
    faceMask = imfill(faceMask, 'holes');
    faceMask = ExtractNLargestBlobs(faceMask, 1);

    faceMask = imdilate(faceMask, strel('disk', 8));
    faceMask = imerode(faceMask, strel('disk', 10));
%     figure; imshow(faceMask); title('fill largest dilate -> erotion (closing) faceMask'); pause; 
    
    faceMask = shrinkFaceMask(faceMask);
%     figure; imshow(faceMask); title('shrink faceMask'); pause; 
    faceMaskRep = repmat(faceMask, [1,1,3]);
    face = rgbImage.*uint8(faceMaskRep);

    grayFace = extractGrayFace(Y, face, rgbImage, faceMaskRep);

%     figure; imshow(grayFace); title('grayFace'); pause;
    grayFace = im2bw(grayFace, 0.47);
%     figure; imshow(grayFace); title('threshold grayFace'); pause;
    


    filteredFaceMask = imfilter(grayFace, fspecial('laplacian'));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause
    filteredEyeCandidates1 = filteredFaceMask;
    filteredFaceMask = bwareaopen(filteredFaceMask, 40);
    originalFilteredFaceMask = filteredFaceMask;
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
   
    [labeledImage, numberOfBlobs] = bwlabel(filteredFaceMask);
    convexProperties = regionprops(labeledImage, 'ConvexArea');
    convexAreas = [convexProperties.ConvexArea];
    [sortedConvexAreas, sortedIndices] = sort(convexAreas, 2, 'descend');
    
    filteredFaceMask = ismember(labeledImage, sortedIndices(1));
%     figure; imshow(filteredFaceMask); title('largest convex area filteredFaceMask'); pause;
    filteredFaceMaskFilled = imfill(filteredFaceMask, 'holes');
%     figure; imshow(filteredFaceMaskFilled); title('filteredFaceMaskFilled'); pause;
    filteredFaceMaskFilledEroded = imerode(filteredFaceMaskFilled, strel('disk', 1));
%     figure; imshow(filteredFaceMaskFilledEroded); title('filteredFaceMaskFilledEroded'); pause;
    filteredFaceMaskBorder = xor(filteredFaceMaskFilled, filteredFaceMaskFilledEroded);
%     figure; imshow(filteredFaceMaskBorder); title('xor(filteredFaceMaskFilled, filteredFaceMaskFilledEroded)'); pause;
    
    originalFilteredFaceMask = imdilate(originalFilteredFaceMask, strel('disk', 1));
    originalFilteredFaceMask = imerode(originalFilteredFaceMask, strel('disk', 1));
    filteredEyeCandidates2 = originalFilteredFaceMask - and(originalFilteredFaceMask, filteredFaceMaskBorder); 
    filteredEyeCandidates2 = imfill(filteredEyeCandidates2, 'holes');
%     figure; imshow(filteredEyeCandidates2); title('filteredEyeCandidates2'); pause;
    
%     filla = imfill(filteredFaceMaskBorder, 'holes');
%     filla = imerode(filla, strel('disk', 5));
%     filteredFaceMaskCopy3 = and(filteredFaceMaskCopy2, filla);
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
%     figure; imshow(filteredFaceMaskCopy2); title('filteredFaceMask'); pause;
    
%     filteredFaceMaskEdyes = filteredFaceMaskFilled;
%     figure; imshow(filteredFaceMaskEyes); title('filteredFaceMaskEyes'); pause;
    
    filteredEyeCandidates1 = filteredEyeCandidates1 - filteredFaceMask;
    filteredEyeCandidates1 = imdilate(filteredEyeCandidates1, strel('disk', 5));
    filteredEyeCandidates1 = imfill(filteredEyeCandidates1, 'holes');
%     figure; imshow(filteredEyeCandidates1); title('filteredEyeCandidates1'); pause;
    
    filteredFaceMaskEyesReal = and(filteredFaceMaskFilled, filteredEyeCandidates1);
    filteredFaceMaskEyesReal = and(filteredFaceMaskEyesReal, filteredEyeCandidates2);
%     figure; imshow(filteredFaceMaskEyesReal); title('filteredFaceMaskEyesReal'); pause;

    
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk', 30));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMask2 = imfill(filteredFaceMask, 'holes');
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMask = ~xor(~filteredFaceMask, filteredFaceMask2);
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMask = ExtractNLargestBlobs(filteredFaceMask, 1);
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    
%     faceMask = and(filteredFaceMask, faceMask);
%     faceMaskRep = repmat(faceMask, [1,1,3]);
%     face = rgbImage2.*uint8(faceMaskRep);
%     figure; imshow(face); title('face'); pause;
end


