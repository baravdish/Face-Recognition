function [normalizedImg] = otsuNormalize(img)

    [bins, binIndexVec] = imhist(img);

    nPixels = length(img(:));
    probabilityBins = bins/nPixels;
    nBins = length(binIndexVec);
    
    WbVec = zeros(nBins, 1);
    WfVec = zeros(nBins, 1);
    sumBVec = zeros(nBins, 1);
    
    for k = 1:nBins
        WbVec(k) = sum(bins(1:k));
        sumBVec(k) = binIndexVec(1:k)'*bins(1:k);
    end

    WfVec = nPixels - WbVec;

    sumTerm = sum(binIndexVec'*bins);
    meanBack = sumBVec./WbVec;
    meanForground = (sumTerm - sumBVec)./WfVec;
    inBetween = WbVec.*WfVec.*(meanBack - meanForground).^2;

    [~, threshold] = max(inBetween);

    normalizedImg = img > threshold;

end