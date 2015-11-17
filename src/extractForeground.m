function foregroundMask = extractForeground(grayImage, overSaturatedMask, rgbImage)
    
    grayImageWithoutOverSaturation = grayImage.*uint8(~overSaturatedMask);
    
    adjustedGrayImageWithoutOversaturation = imadjust(grayImageWithoutOverSaturation);
    
    lowFrequencyParts = medfilt2(adjustedGrayImageWithoutOversaturation);
    highFrequencyParts = imfilter(lowFrequencyParts, fspecial('log'));
    
%     imshow(highFrequencyParts); title('highFrequencyParts'); pause;
    
% %     Try with otsu threshold
%     otsuLevel = graythresh(rgbImage);
%     bw = im2bw(rgbImage, otsuLevel);
%     bw = bwareaopen(bw, 50);
%     bw = ~bw;
%     bw = imfill(bw, 'holes');
% 
%     surface1 = highFrequencyParts.*uint8(bw);
%     surface2 = highFrequencyParts.*uint8(~bw);
%     
%     surface1HighFrequencyParts = highFrequencyParts.*uint8(surface1);
%     foregroundPercentage = (100 * sum(surface1HighFrequencyParts(:))) / (length(nonzeros(surface1(:)))*255);
%     
%     surface2HighFrequencyParts = highFrequencyParts.*uint8(surface2);
%     backgroundPercentage = (100 * sum(surface2HighFrequencyParts(:))) / (length(nonzeros(surface2(:)))*255);
%     
%     foreground = bw;
%       
%     if backgroundPercentage > foregroundPercentage
%        foreground = ~bw; 
%        temp = foregroundPercentage;
%        foregroundPercentage = backgroundPercentage;
%        backgroundPercentage = temp;
%     end
%     
%     if backgroundPercentage < 30 && (foregroundPercentage / backgroundPercentage) > 2
%         % The large background found by otsu is homogenous
%         foregroundMask = imdilate(foreground, strel('disk', 10));
%         foregroundMask = imfill(foregroundMask, 'holes');
%         foregroundMask = imerode(foregroundMask, strel('disk', 10));
% %         figure; imshow(foregroundMask); title('foregroundMask by otsu'); pause;
%         return;
%     end
    
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