function binaryImage = extractNLargestBlobs(binaryImage, numberToExtract)
    [labeledImage, numberOfBlobs] = bwlabel(binaryImage);
	blobMeasurements = regionprops(labeledImage, 'area');

	allAreas = [blobMeasurements.Area];
	if numberToExtract > 0
		[sortedAreas, sortIndexes] = sort(allAreas, 'descend');
	elseif numberToExtract < 0
		[sortedAreas, sortIndexes] = sort(allAreas, 'ascend');
		numberToExtract = -numberToExtract;
	else

		binaryImage = false(size(binaryImage));
		return;
    end
    
	biggestBlob = ismember(labeledImage, sortIndexes(1:numberToExtract));
	
	binaryImage = biggestBlob > 0;

end