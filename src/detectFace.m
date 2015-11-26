function output = detectFace(rgbImage)
%     figure; imshow(rgbImage); title('rgbImage'); pause;   
    
    rgbImage = padarray(rgbImage,[2 2],'both');
    
    grayImage = rgb2gray(rgbImage);
%     figure; imshow(grayImage); title('grayImage'); pause;
    
    [estimatedSkinMask, Y, Cb, Cr] = extractSkinColorFromYCbCr(rgbImage);
%     figure; imshow(estimatedSkinMask); title('estimatedSkinMask'); pause;   
    
    overSaturatedMask = extractOverSaturated(grayImage);
%     figure; imshow(overSaturatedMask); title('overSaturatedMask'); pause;   
    
    foregroundMask = extractForeground(grayImage, overSaturatedMask, rgbImage);
%     figure; imshow(foregroundMask); title('foregroundMask'); pause; 
    
    faceMask = and(foregroundMask, estimatedSkinMask);
%     figure; imshow(faceMask); title('faceMask'); pause;
    
    nonSkinMask1 = Cb > Cr + 10;
%     figure; imshow(nonSkinMask1); title('nonSkinMask1'); pause;
    
    hsvImage = rgb2hsv(rgbImage);
    H = hsvImage(:,:,1);
    S = hsvImage(:,:,2);
    V = hsvImage(:,:,3);
    
    nonSkinMask2 = V < S - 0.2;
%     figure; imshow(nonSkinMask2); title('nonSkinMask2'); pause;

%     hairMask = nonSkinMask2;
%     figure; imshow(hairMask); title('hairMask'); pause;
%     hairMask = imdilate(hairMask, strel('disk', 2));
%     hairMask = ExtractNLargestBlobs(hairMask, 1);
%     hairMask = imerode(hairMask, strel('disk', 2));
%     
% %     nonSkinMask3 = bwareaopen(nonSkinMask3, 1000);
%     figure; imshow(hairMask); title('hairMask'); pause;    
    
%     nonSkinMask3 = V < S;
%     nonSkinMask3 = imdilate(nonSkinMask3, strel('disk', 2));
%     nonSkinMask3 = imerode(nonSkinMask3, strel('disk', 2));
%     nonSkinMask3 = ~nonSkinMask3;
%     figure; imshow(nonSkinMask3); title('nonSkinMask3'); pause;
%     nonSkinMask4 = imfill(nonSkinMask3, 'holes');
%     figure; imshow(nonSkinMask4); title('nonSkinMask4'); pause;
%     nonSkinMask3 = xor(nonSkinMask3, nonSkinMask4);
%     figure; imshow(nonSkinMask3); title('nonSkinMask3'); pause;
%     nonSkinMask3 = imfill(nonSkinMask3, 'holes');
%     nonSkinMask3 = imdilate(nonSkinMask3, strel('disk', 8));
%     figure; imshow(nonSkinMask3); title('nonSkinMask3'); pause;
    
    noFaceMask =  or(nonSkinMask1, nonSkinMask2);
    noFaceMask = imerode(noFaceMask, strel('disk', 5));
    noFaceMask = bwareaopen(noFaceMask, 1000);
    noFaceMask = imdilate(noFaceMask, strel('disk', 5));
%     figure; imshow(noFaceMask); title('noFaceMask'); pause;
    
    noFaceMaskRep = repmat(noFaceMask, [1,1,3]);
    noFace = rgbImage.*uint8(noFaceMaskRep);
%     figure; imshow(noFace); title('noFace'); pause;
    
    faceMask = and(faceMask, ~noFaceMask);
%     figure; imshow(faceMask); title('faceMask'); pause;
    
    faceMask = imfill(faceMask, 'holes');
    faceMask = ExtractNLargestBlobs(faceMask, 1);
%     figure; imshow(faceMask); title('faceMask'); pause;

    faceMask = imdilate(faceMask, strel('disk', 8));
    faceMask = imerode(faceMask, strel('disk', 10));
    
%     figure; imshow(faceMask); title('faceMask'); pause;
    
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
    
    faceMaskRep = repmat(faceMask, [1,1,3]);
    face = rgbImage.*uint8(faceMaskRep);
%     figure; imshow(face); title('face'); pause;
    
    filteredFaceMask = imfilter(im2bw(imadjust(rgb2gray(face)), 0.47), fspecial('laplacian'));
    filteredFaceMaskCopy = filteredFaceMask;
    
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    
    filteredFaceMask = bwareaopen(filteredFaceMask, 40);
    originalFilteredFaceMask = filteredFaceMask;
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
   
    [labeledImage, numberOfBlobs] = bwlabel(filteredFaceMask);
    convexProperties = regionprops(labeledImage, 'ConvexArea');
    convexAreas = [convexProperties.ConvexArea];
    [sortedConvexAreas, sortedIndices] = sort(convexAreas, 2, 'descend');
    
    filteredFaceMask = ismember(labeledImage, sortedIndices(1));
    
    originalFilteredFaceMask = imdilate(originalFilteredFaceMask, strel('disk', 1));
    originalFilteredFaceMask = imerode(originalFilteredFaceMask, strel('disk', 1));
%     figure; imshow(originalFilteredFaceMask); title('originalFilteredFaceMask'); pause;
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMaskCopy2 = originalFilteredFaceMask - and(originalFilteredFaceMask, filteredFaceMask); 
%     figure; imshow(filteredFaceMaskCopy2); title('filteredFaceMaskCopy2'); pause;
    filteredFaceMaskCopy2 = imfill(filteredFaceMaskCopy2, 'holes');
%     figure; imshow(filteredFaceMaskCopy2); title('filteredFaceMaskCopy2'); pause;
    
    filteredFaceMaskEyes = filteredFaceMask;
    filteredFaceMaskEyes = imfill(filteredFaceMaskEyes, 'holes');
%     figure; imshow(filteredFaceMaskEyes); title('filteredFaceMaskEyes'); pause;
    filteredFaceMaskCopy = filteredFaceMaskCopy - filteredFaceMask;
    filteredFaceMaskCopy = imfill(filteredFaceMaskCopy, 'holes');
    filteredFaceMaskCopy = imdilate(filteredFaceMaskCopy, strel('disk', 5));
    filteredFaceMaskCopy = imfill(filteredFaceMaskCopy, 'holes');
%     
%     figure; imshow(filteredFaceMaskCopy); title('filteredFaceMaskCopy'); pause;
%     filteredFaceMaskCopy2 = originalFilteredFaceMask - and(originalFilteredFaceMask, filteredFaceMaskCopy); 
%     figure; imshow(filteredFaceMaskCopy2); title('filteredFaceMaskCopy2'); pause;
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk', 30));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;

    filteredFaceMask2 = filteredFaceMask;
    filteredFaceMask2 = imfill(filteredFaceMask2, 'holes');
%     figure; imshow(filteredFaceMask2); title('filteredFaceMask2'); pause;
    filteredFaceMask = ~xor(~filteredFaceMask, filteredFaceMask2);
%     figure; imshow(filteredFaceMask); title('filteredFaceMask5'); pause;

    filteredFaceMask = ExtractNLargestBlobs(filteredFaceMask, 1);
%     figure; imshow(filteredFaceMask); title('filteredFaceMask7'); pause;
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask9'); pause;
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask10'); pause;
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask11'); pause;
    
    faceMask = and(filteredFaceMask, faceMask);
    faceMaskRep = repmat(faceMask, [1,1,3]);
    face = rgbImage.*uint8(faceMaskRep);
%     figure; imshow(face); title('face'); pause;
    
    
    nonSkinMask3 = meanColorize(rgbImage, faceMask, 0.2);
    nonSkinMask3 = imerode(nonSkinMask3, strel('disk', 2));
    nonSkinMask3 = imfill(nonSkinMask3, 'holes');

















    
    
    

    
    
    
    
    
  
%     Cr = im2double(Cr);
%     Cb = im2double(Cb);
    
%     max(Cr(:))
%     min(Cr(:))
%     max(Cb(:))
%     min(Cb(:))
    
    Cb_orig = Cb;
    Cr_orig = Cr;
    Y = uint16(Y);
    Cr = uint16(Cb);
    Cb = uint16(Cb);
    
    CbSquared = Cb.^2;
        
    CrInvertedSquared = (1-Cr).^2;
    CbDividedByCr = Cb./Cr;
    
    eyeMapChroma = (1/3)*(CbSquared + ...
                          CrInvertedSquared + ...
                          CbDividedByCr);
%     max(eyeMapChroma(:))
%     min(eyeMapChroma(:))
                      
                      
    eyeMapLuma = uint16(zeros(size(grayImage)));
%     eyeMapLuma = im2double(grayImage);
%     grayImage = imadjust(grayImage);
    for n=1 : 10
        kernel = strel('disk', n);
        dilation = imdilate(Y, kernel);
        erosion = imerode(Y, kernel);
        eyeMapLuma = eyeMapLuma + (dilation./(n + erosion));
%         figure; imshow(eyeMapLuma); title('eyeMapLuma'); pause;
    end
%     max(eyeMapLuma(:))
%     min(eyeMapLuma(:))
    
%     eyeMapLuma = imadjust(eyeMapLuma);
    
    grayFace = im2double(imadjust(rgb2gray(face)));
    filtered = im2double(zeros(size(grayFace)));
    for k = 1:10
        h = (ones(k, k) - fspecial('gaussian', k, 0.1))./k^2;
        filtered = filtered + im2double(imfilter(grayFace, h)) / k;
        filtered = imadjust(filtered);
    end
    filtered = im2double(ones(size(filtered))) - filtered;
    filtered = filtered.*faceMask;
%     eyeMapLuma = imdilate(eyeMapLuma, strel('disk', 1));
    
    eyeMapChroma = histeq(eyeMapChroma);
    eyeMapChroma = im2double(eyeMapChroma);
%     eyeMapChroma = 255 * ((eyeMapChroma - min(eyeMapChroma(:))) ...
%                  ./ (max(eyeMapChroma(:)) - min(eyeMapChroma(:))));

    
%     figure; imshow(eyeMapChroma, []); title('eyeMapChroma'); pause;
    
    eyeMapLuma = im2double(eyeMapLuma);
%     eyeMapLuma = 255 * ((eyeMapLuma - min(eyeMapLuma(:))) ...
%                  ./ (max(eyeMapLuma(:)) - min(eyeMapLuma(:))));

%     max(eyeMapLuma(:))
%     min(eyeMapLuma(:))
    
%     figure; imshow(eyeMapLuma, []); title('eyeMapLuma'); pause;
    
    eyeMap = eyeMapChroma .* eyeMapLuma;
    
    eyeMap = im2double(eyeMap);
    eyeMap = 255 * ((eyeMap - min(eyeMap(:))) ...
                 ./ (max(eyeMap(:)) - min(eyeMap(:))));
    eyeMap = im2double(eyeMap);
             
    originalEyeMap = eyeMap;
    
%     figure; imshow(eyeMap, []); title('originalEyeMap'); pause;
    
    eyeMap = eyeMap.*filteredFaceMaskCopy;
%     figure; imshow(filteredFaceMaskCopy); title('filteredFaceMaskCopy'); pause;
%     figure; imshow(eyeMap); title('eyeMap.*filteredFaceMaskCopy'); pause;
  
    eyeMap = eyeMap.*filteredFaceMaskCopy2;
%     figure; imshow(filteredFaceMaskCopy2); title('filteredFaceMaskCopy2'); pause;
%     figure; imshow(eyeMap); title('eyeMap.*filteredFaceMaskCopy2'); pause;
    
    eyeMap = eyeMap.*im2double(faceMask);
%     figure; imshow(eyeMap); title('eyeMap.*im2double(faceMask)'); pause;
    



    
    
    
   
%     eyeMap = eyeMap - im2double(mouthMap);

    
    
%     figure; imshow(eyeMap); title('eyeMap11'); pause;
    eyeMap = eyeMap - and(overSaturatedMask, eyeMap);
%     figure; imshow(eyeMap); title('eyeMap2'); pause;
    eyeMap = imfill(eyeMap, 'holes');
    eyeMap = imdilate(eyeMap, strel('disk', 4));
%     figure; imshow(eyeMap); title('before final Eyemap'); pause;
    
%     eyeMap = ExtractNLargestBlobs(eyeMap, 2);
% %     figure; imshow(eyeMap); title('final Eyemap'); pause;
% %     figure; imshow(eyeMap); title('final Eyemap'); pause;
%     eyeMap = imerode(eyeMap, strel('disk', 3));
% %     figure; imshow(eyeMap); title('final Eyemap'); pause;
%     eyeMap = imdilate(eyeMap, strel('disk', 14));
% %     figure; imshow(eyeMap); title('final Eyemap'); pause;
%     eyeMap = imfill(eyeMap, 'holes');
% %     figure; imshow(eyeMap); title('final Eyemap'); pause;
%     eyeMap = ExtractNLargestBlobs(eyeMap, 2);
% %     figure; imshow(eyeMap); title('final Eyemap'); pause;
    
    eyeMapRep = repmat(eyeMap, [1,1,3]);
    I = face.*uint8(eyeMapRep);
    
%     figure; imshow(estimatedSkinMask); title('estimatedSkinMask'); pause;
    estimatedSkinMask = imerode(estimatedSkinMask, strel('disk', 5));
%     figure; imshow(estimatedSkinMask); title('estimatedSkinMask'); pause;
%     figure; imshow(eyeMap); title('eyeMap'); pause;
    estimatedSkinMask = imerode(estimatedSkinMask, strel('disk', 5));
    eyeMap = eyeMap .* ~estimatedSkinMask;
    
%     figure; imshow(eyeMap); title('eyeMap'); pause;

    eyeMap = imfill(eyeMap, 'holes');
    eyeMap = eyeMap > 0; 
    eyeMap = imdilate(eyeMap, strel('disk', 10));

    eyeMapRep = repmat(eyeMap, [1,1,3]);
    I = rgbImage.*uint8(eyeMapRep);
%     I = rgbImage;
% figure; imshow(originalEyeMap); title('originalEyeMap'); pause;
    IGradient = imgradient(originalEyeMap);
    IGradient = originalEyeMap;
%      figure; imshow(IGradient, []); title('IGradient'); pause;
%     IGradient = imfilter(IGradient, fspecial('log'));
%     IGradient = edge(IGradient,'Canny');
%      figure; imshow(IGradient, []); title('IGradient'); pause;
     
    IGradient = IGradient.*eyeMap;
%     figure; imshow(IGradient); title('IGradient'); pause;
    IGradient = IGradient.*filteredFaceMaskEyes;
%     figure; imshow(IGradient); title('IGradient'); pause;
%     IGradient = IGradient.*skinH;
%     IGradient = IGradient.*filtered;
%     faceMask = imerode(faceMask, strel('disk',4));
    IGradient = IGradient.*faceMask;
  
    
%     figure; imshow(IGradient, []); title('IGradient'); pause;

    
%     figure; imshow(nonSkinMask3); title('nonSkinMask3'); pause;
    IGradient = IGradient.*nonSkinMask3;
%     figure; imshow(IGradient, []); title('IGradient'); pause;
    
%     nonSkinMask1 = Cb >= Cr;
%     nonSkinMask1 = imdilate(nonSkinMask1, strel('disk', 10));
%     IGradient = IGradient.*nonSkinMask1;
%     IGradient = IGradient.*meanColorizedFace;
%     figure; imshow(nonSkinMask1); title('nonSkinMask1'); pause;
%     figure; imshow(I); title('I dark'); pause;
    
    IBright = IGradient;
    IDark = IGradient;
    
%      figure; imshow(IBright, []); title('IBright'); pause;

    [centersBright, radiiBright] = imfindcircles(IBright, [4 20], ...
                 'ObjectPolarity', 'bright', 'Method', 'TwoStage', ...
             'sensitivity', 0.99);
             
%     figure; imshow(I); title('I'); pause;
%     viscircles(centersDark(1:2,:), radiiDark(1:2),'EdgeColor','r');
%     viscircles(centersBright(1:2,:), radiiBright(1:2),'EdgeColor', 'b');
    
%     viscircles(centersDark(1:10,:), radiiDark(1:10), 'EdgeColor', 'r');
%     viscircles(centersDark(1:2,:), radiiDark(1:2),'EdgeColor', 'g');
%     viscircles(centersBright(1:10,:), radiiBright(1:10),'EdgeColor', 'b');

% pause;

%     
    centersDark = centersBright;
    radiiDark = radiiBright;
    
    eyeMapSize = size(I);
    imageSizeX = eyeMapSize(1,1,1);
    imageSizeY = eyeMapSize(1,2,1);
    [cols rows] = meshgrid(1:imageSizeY, 1:imageSizeX);

    centerX1 = centersDark(1,1,1);
    centerY1 = centersDark(1,2,1);
    radius1 = radiiDark(1,1);
    irisMask = ((rows - centerY1).^2 + (cols - centerX1).^2) <= radius1.^2;
    
    centerX2 = centersDark(2,1,1);
    centerY2 = centersDark(2,2,1);
    radius2 = radiiDark(2,1);
    
    for n=2 : size(centersDark, 1)
        centerX2 = centersDark(n,1,1);
        centerY2 = centersDark(n,2,1) ; 
        radius2 = radiiDark(n,1);
        
        distance = sqrt((centerX2-centerX1)^2 + (centerY2-centerY1)^2);
        
        if distance > (radius1 + radius2)*4
            n = size(centersDark, 1);
            break
        end
    end
         

    irisMask = irisMask | ((rows - centerY2).^2 + (cols - centerX2).^2) <= radius2.^2;
    
    irisMaskRep = repmat(irisMask, [1,1,3]);
    irisMap = face.*uint8(irisMaskRep);
    
%     figure; imshow(irisMap); title('irisMap'); pause;
    
    leftEyeCenterX = centerX1;
    leftEyeCenterY = centerY1;
    leftRadius = radius1;
    
    rightEyeCenterX = centerX2;
    rightEyeCenterY = centerY2;
    rightRadius = radius2;
    if centerX1 > centerX2 
        leftEyeCenterX = centerX2;
        leftEyeCenterY = centerY2;
        leftRadius = radius2;
        
        rightEyeCenterX = centerX1;
        rightEyeCenterY = centerY1;  
        rightRadius = radius1;
    end
    
    
    
    eyeDistance = rightEyeCenterX - leftEyeCenterX;
    measureDistance = eyeDistance * 0.9;
    eyeY = (leftEyeCenterY + rightEyeCenterY) / 2;

    eyeRegionMask = zeros(size(face(:,:,1)));
%     mask(eyeY-measureDistance : eyeY+measureDistance, :) = 1;
    eyeRegionMask(1: eyeY+measureDistance, :) = 1;
%     figure; imshow(mask); title('mask'); pause;
    
    

%     Cb = im2double(Cb);
%     Cr = im2double(Cr);
    
%     n = length(face(:,:,1));
    n = length(nonzeros(faceMask));
    
    Cb = double(Cb_orig);
    Cr = double(Cr_orig);
    
    CrSquared = Cr.^2;
    numerator = (1 / n) * sum(CrSquared(:));
    
    CrDividedByCb = Cr./Cb;
    denumerator = (1 / n) * sum(CrDividedByCb(:));
    
    nn = 0.95 * (numerator / denumerator);
    
    mouthMap = CrSquared.*(CrSquared - nn*CrDividedByCb).^2;
    
%     figure; imshow(mouthMap, []); title('mouthMap'); pause;
    
    mouthMap = im2double(mouthMap);
    mouthMap = 255 * ((mouthMap - min(mouthMap(:))) ...
                 ./ (max(mouthMap(:)) - min(mouthMap(:))));
    mouthMap = uint8(mouthMap);
    
    originalMouthMap = mouthMap;
    
%     figure; imshow(mouthMap, []); title('mouthMap'); pause;
   
    
    max(mouthMap(:))
    min(mouthMap(:))
    

%     figure; imshow(mouthMap, []); title('mouthMap'); pause;
    mouthMap = mouthMap .* uint8(faceMask);
%     figure; imshow(nonSkinMask3, []); title('nonSkinMask3'); pause;
%     figure; imshow(estimatedSkinMask, []); title('estimatedSkinMask'); pause;
%     estimatedSkinMask2 = ExtractNLargestBlobs(estimatedSkinMask, 1);
%     figure; imshow(estimatedSkinMask2, []); title('estimatedSkinMask2'); pause;
%     estimatedSkinMask2 = imerode(estimatedSkinMask2, strel('disk', 3));
% %     estimatedSkinMask2 = imdilate(estimatedSkinMask2, strel('disk', 3));
%     figure; imshow(estimatedSkinMask2, []); title('estimatedSkinMask2'); pause;
%     estimatedSkinMask2 = imfill(estimatedSkinMask2, 'holes');
%     figure; imshow(estimatedSkinMask2, []); title('estimatedSkinMask2'); pause;
%     estimatedSkinMask2 = imerode(estimatedSkinMask2, strel('disk', 2));
%     mouthMap = mouthMap .* uint8(estimatedSkinMask2);
%     figure; imshow(mouthMap, []); title('mouthMap'); pause;
    
%     figure; imshow(face, []); title('face'); pause;
%     figure; imshow(mouthMap, []); title('mouthMap'); pause;
%     mouthMap

    mouthMap = mouthMap > 0.5 * max(mouthMap(~eyeRegionMask));

%     meana = mean2(mouthMap(faceMask));
%     maxa = max(mouthMap(:));
%     mouthMap = mouthMap > meana + 0.15*(maxa - meana);

%     mouthMap = im2bw(mouthMap, graythresh(mouthMap));
    
%     figure; imshow(mouthMap, []); title('mouthMap'); pause;
    mouthMap = imdilate(mouthMap, strel('disk', 4));
    mouthMap = imfill(mouthMap, 'holes');
    mouthMap = imerode(mouthMap, strel('disk', 4));
%     figure; imshow(mouthMap, []); title('mouthMap'); pause;
%     mouthMap = rgb2gray(face) .* uint8(mouthMap);
%     mouthMap = originalMouthMap .* uint8(mouthMap);
%     mouthMapRep = repmat(mouthMap, [1,1,3]);
%     mouthMap = face.*uint8(mouthMapRep);
    
%     figure; imshow(mouthMap, []); title('mouthMap'); pause;
%     mouthMap = or(imfilter(mouthMap, [-1 0 1]'), imfilter(mouthMap, [1 0 -1]'));
%     figure; imshow(mouthMap, []); title('mouthMap'); pause;
%     mouthMap = imdilate(mouthMap, strel('disk', 2));
%     mouthMap = imerode(mouthMap, strel('disk', 2));
%     figure; imshow(mouthMap); title('mouthMap'); pause;
%     mouthMap = imfill(mouthMap, 'holes');
%     figure; imshow(mouthMap); title('mouthMap'); pause;   
    



    

    BW = and(mouthMap, ~eyeRegionMask);
%     BW = mouthMap > 0.8;
%     figure; imshow(BW); title('BW'); pause;
    
    BW = ExtractNLargestBlobs(BW, 1);
%     figure; imshow(BW); title('Binary Image'); pause;
    s = regionprops(BW, mouthMap, {'Centroid','WeightedCentroid', 'Area'});
    
%     figure;
%     imshow(mouthMap)
%     title('Weighted (red) and Unweighted (blue) Centroids');
%     hold on
%     numObj = numel(s);
%     for k = 1 : numObj
%         plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
%         plot(s(k).Centroid(1), s(k).Centroid(2), 'bo');
%     end
%     hold off
%     pause;
    mouthX = s(1).Centroid(1);
    mouthY = s(1).Centroid(2);

    eyeMouthDistance = mouthY - eyeY;
    
    offsetX = eyeDistance * 0.3
    offsetY = eyeMouthDistance * 0.2;
    
    output = rgbImage(round(eyeY-offsetY) : round(mouthY+offsetY), ...
                      round(leftEyeCenterX-offsetX) :  round(rightEyeCenterX+offsetX), :);
    figure; imshow(output); title('output');
    
    
    
%     figure; imshow(face); title('face before rotated'); 
    viscircles([leftEyeCenterX, leftEyeCenterY], leftRadius,'EdgeColor','r');
    viscircles([rightEyeCenterX, rightEyeCenterY], rightRadius,'EdgeColor', 'g');
    viscircles([mouthX, mouthY], rightRadius,'EdgeColor', 'b');
%     pause;
    
    referenceDirection = [0 1];
    eyeDirection = [leftEyeCenterX-rightEyeCenterX, leftEyeCenterY-rightEyeCenterY];
    eyeDirection = eyeDirection / norm(eyeDirection);
    
    angle = acos(dot(referenceDirection, eyeDirection)) * 180 / pi;
    angle = 90 - angle;

%     figure; imshow(face); title('face before rotated'); pause;
    
    marker=zeros(face(1:2));
    
    face = imrotate(face, -angle);
    
    marker(round(leftEyeCenterX), round(leftEyeCenterY))=1;
    marker_rot = imrotate(marker, -angle);
    [leftEyeCenterX, leftEyeCenterY]=find(marker_rot)
    
    s=size(face);
    marker=zeros(s(1:2));
    marker(round(rightEyeCenterX), round(rightEyeCenterY))=1;
    marker_rot = imrotate(marker, -angle);
    [rightEyeCenterX, rightEyeCenterY]=find(marker_rot)

    
%     [rightEyeCenterX, rightEyeCenterY] = imrotate([rightEyeCenterX, rightEyeCenterY], -angle);
    
%     figure; imshow(face); title('face after rotated'); pause;
%     viscircles([leftEyeCenterX, leftEyeCenterY], leftRadius,'EdgeColor','r');
%     viscircles([rightEyeCenterX, rightEyeCenterY], rightRadius,'EdgeColor', 'g');
%     pause;
    
    output = rgbImage;
end
