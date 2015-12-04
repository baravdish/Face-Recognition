function foregroundMask = extractForeground(grayImage, overSaturatedMask, rgbImage)
    
    grayImageWithoutOverSaturation = grayImage.*uint8(~overSaturatedMask);
    
    adjustedGrayImageWithoutOversaturation = imadjust(grayImageWithoutOverSaturation);
    
    lowFrequencyParts = medfilt2(adjustedGrayImageWithoutOversaturation);
    highFrequencyParts = imfilter(lowFrequencyParts, fspecial('log'));
        
    edges = edge(adjustedGrayImageWithoutOversaturation, 'canny');

    blobs = ~imdilate(edges, strel('disk', 2));
    blobs = imerode(blobs, strel('disk', 2));

%     figure; imshow(blobs); title('All background blobs'); pause;
    
    largestBlobSize = size(grayImage, 1) * size(grayImage, 2) * 0.05;
    blobs = bwareaopen(blobs, floor(largestBlobSize));
    
%     figure; imshow(blobs); title('Large background blobs'); pause;

    [labeledImage, numberOfBlobs] = bwlabel(blobs);
    blobMeasurements = regionprops(labeledImage, 'area');

    allAreas = [blobMeasurements.Area];
    [sortedAreas, sortIndexes] = sort(allAreas, 'descend');

    background = zeros(size(grayImage));
    
    for n = 1:length(sortIndexes)
        currentBlob = ismember(labeledImage, sortIndexes(n));
%         figure; imshow(currentBlob); title('currentBlob'); pause;
        currentHighFrequencyParts = highFrequencyParts.*uint8(currentBlob);
        percentageOfHighFrequency = (100 * sum(currentHighFrequencyParts(:))) / (length(nonzeros(currentBlob(:)))*255);
%         figure; imshow(currentHighFrequencyParts); title('currentHighFrequencyParts'); pause;
        if percentageOfHighFrequency < 0.8
        	background = background | currentBlob;
        end
    end
    
    foregroundMask = ~background;
    
%     imshow(background); title('Background mask'); pause;
%     foregroundMaskRep = repmat(foregroundMask,[1,1,3]);
%     foreground = rgbImage.*uint8(foregroundMaskRep);
%     imshow(foreground); title('foreground'); pause;

end