function [faceMask, face, ...
          filteredFaceMaskCopy, ...
          filteredFaceMaskCopy2, ...
          filteredFaceMaskEyes] = extractFaceMask(rgbImage, ...
                                                  foregroundMask, ...
                                                  estimatedSkinMask, ...
                                                  nonFaceMask)

    faceMask = foregroundMask & estimatedSkinMask & ~nonFaceMask;
%     figure; imshow(faceMask); title('faceMask'); pause;
   
    faceMask = imfill(faceMask, 'holes');
    faceMask = ExtractNLargestBlobs(faceMask, 1);
    faceMask = imdilate(faceMask, strel('disk', 8));
    faceMask = imerode(faceMask, strel('disk', 10));
%     figure; imshow(faceMask); title('faceMask'); pause; 
    
    faceMask = shrinkFaceMask(faceMask);
%     figure; imshow(faceMask); title('faceMask'); pause; 

    faceMaskRep = repmat(faceMask, [1,1,3]);
    face = rgbImage.*uint8(faceMaskRep);
%     figure; imshow(face); title('face'); pause;

    grayFace = imadjust(rgb2gray(face));
    grayFace = im2bw(grayFace, 0.47);
    
    filteredFaceMask = imfilter(grayFace, fspecial('laplacian'));
    filteredFaceMaskCopy = filteredFaceMask;
    filteredFaceMask = bwareaopen(filteredFaceMask, 40);
    originalFilteredFaceMask = filteredFaceMask;
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
   
    [labeledImage, numberOfBlobs] = bwlabel(filteredFaceMask);
    convexProperties = regionprops(labeledImage, 'ConvexArea');
    convexAreas = [convexProperties.ConvexArea];
    [sortedConvexAreas, sortedIndices] = sort(convexAreas, 2, 'descend');
    
    filteredFaceMask = ismember(labeledImage, sortedIndices(1));
    filteredFaceMaskFilled = imfill(filteredFaceMask, 'holes');
    filteredFaceMaskFilledEroded = imerode(filteredFaceMaskFilled, strel('disk', 1));
    filteredFaceMaskBorder = xor(filteredFaceMaskFilled, filteredFaceMaskFilledEroded);
    
    originalFilteredFaceMask = imdilate(originalFilteredFaceMask, strel('disk', 1));
    originalFilteredFaceMask = imerode(originalFilteredFaceMask, strel('disk', 1));
    filteredFaceMaskCopy2 = originalFilteredFaceMask - and(originalFilteredFaceMask, filteredFaceMaskBorder); 
    filteredFaceMaskCopy2 = imfill(filteredFaceMaskCopy2, 'holes');
%     figure; imshow(filteredFaceMaskCopy2); title('filteredFaceMaskCopy2'); pause;
    
    filteredFaceMaskEyes = filteredFaceMaskFilled;
%     figure; imshow(filteredFaceMaskEyes); title('filteredFaceMaskEyes'); pause;
    
    filteredFaceMaskCopy = filteredFaceMaskCopy - filteredFaceMask;
    filteredFaceMaskCopy = imdilate(filteredFaceMaskCopy, strel('disk', 5));
    filteredFaceMaskCopy = imfill(filteredFaceMaskCopy, 'holes');
%     figure; imshow(filteredFaceMaskCopy); title('filteredFaceMaskCopy'); pause;
    
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk', 30));

    filteredFaceMask2 = imfill(filteredFaceMask, 'holes');
    filteredFaceMask = ~xor(~filteredFaceMask, filteredFaceMask2);

    filteredFaceMask = ExtractNLargestBlobs(filteredFaceMask, 1);
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    
    faceMask = and(filteredFaceMask, faceMask);
    faceMaskRep = repmat(faceMask, [1,1,3]);
    face = rgbImage.*uint8(faceMaskRep);
%     figure; imshow(face); title('face'); pause;
end


function faceMask = shrinkFaceMask(faceMask) 
    faceMaskCopy = faceMask;
    minCol = [];
    maxCol = [];
    minRow = [];
    maxRow = [];
    eraseRadius = 80; % 80
    while( isempty(minCol) || isempty(maxCol) || ... 
           isempty(minRow) || isempty(maxRow) )
        faceMaskCopy = faceMask;
        minCol = [];
        maxCol = [];
        minRow = [];
        maxRow = [];
        faceMaskCopy = imerode(faceMaskCopy, strel('disk', eraseRadius));
        faceMaskCopy = imdilate(faceMaskCopy, strel('disk', eraseRadius));
        [row, col] = find(faceMaskCopy(:,:,1) ~= 0);
        minCol = min(col);
        maxCol = max(col);
        minRow = min(row);
        maxRow = max(row);
        eraseRadius = max(eraseRadius - 10, 1);
%         imshow(faceMaskCopy); title('faceMaskCopy'); pause;
    end
    faceMask = faceMaskCopy;
end