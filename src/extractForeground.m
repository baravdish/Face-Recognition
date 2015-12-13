function foregroundMask = extractForeground(grayImage, overSaturatedMask)
    
    grayImageWithoutOverSaturation = grayImage.*uint8(~overSaturatedMask);
    
    adjustedGrayImageWithoutOversaturation = imadjust(grayImageWithoutOverSaturation);
    
    lowFrequencyParts = medfilt2(adjustedGrayImageWithoutOversaturation);
    highFrequencyParts = imfilter(lowFrequencyParts, fspecial('log'));
        
    edges = edge(adjustedGrayImageWithoutOversaturation, 'canny');

    blobs = ~imdilate(edges, strel('disk', 2));
    blobs = imerode(blobs, strel('disk', 2));

    largestBlobSize = size(grayImage, 1) * size(grayImage, 2) * 0.05;
    blobs = bwareaopen(blobs, floor(largestBlobSize));
    
    [labeledImage, numberOfBlobs] = bwlabel(blobs);
    blobMeasurements = regionprops(labeledImage, 'area');

    allAreas = [blobMeasurements.Area];
    [sortedAreas, sortIndexes] = sort(allAreas, 'descend');

    background = zeros(size(grayImage));
    
    for n = 1:length(sortIndexes)
        currentBlob = ismember(labeledImage, sortIndexes(n));
        currentHighFrequencyParts = highFrequencyParts.*uint8(currentBlob);
        percentageOfHighFrequency = (100 * sum(currentHighFrequencyParts(:))) / (length(nonzeros(currentBlob(:)))*255);
        if percentageOfHighFrequency < 0.8
        	background = background | currentBlob;
        end
    end
    
    foregroundMask = ~background;

end