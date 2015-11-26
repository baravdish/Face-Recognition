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
    
    noFaceMask =  nonSkinMask1 | nonSkinMask2;
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
    eraseRadius = 80;
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
    
    filteredFaceMask = bwareaopen(filteredFaceMask, 50);
    originalFilteredFaceMask = filteredFaceMask;
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
   
    [labeledImage, numberOfBlobs] = bwlabel(filteredFaceMask);
    convexProperties = regionprops(labeledImage, 'ConvexArea');
    convexAreas = [convexProperties.ConvexArea];
    [sortedConvexAreas, sortedIndices] = sort(convexAreas, 2, 'descend');
    
    filteredFaceMask = ismember(labeledImage, sortedIndices(1));
    
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMaskCopy2 = originalFilteredFaceMask - and(originalFilteredFaceMask, filteredFaceMask); 
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
%     figure; imshow(filteredFaceMask); title('filteredFaceMask4'); pause;

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
    
    
    
    
    
    
    
%     
%     figure; imshow(faceMask); title('faceMask'); pause;
%     [labeledFaceMask, numberOfBlobsInFaceMask] = bwlabel(faceMask);
%     [props] = regionprops(labeledFaceMask, 'Orientation');
%     orientation = props.Orientation;
% 
%     orientation = 90-orientation
%     if orientation > 90
%         orientation = orientation - 180
%     end
%     
%     figure; imshow(face); title('face before rotated'); pause;
%     face = imrotate(face, orientation);
%     figure; imshow(face); title('face after rotated'); pause;
%     
%     output = rgbImage;
%     return
%     
    
    
    
    
    
    
    
    
    

%     figure; imshow(face); title('face'); pause;
  
%     Y = im2double(Y);
    Cr = im2double(Cr);
    Cb = im2double(Cb);
    
    CbSquared = Cb.^2;
    CrInvertedSquared = (1-Cr).^2;
%     CbSquaredDiviedByCr = CbSquared./Cr;
    CbDividedByCr = Cb./Cr;
    eyeMapChroma = (1/3)*(CbSquared + ...
                          CrInvertedSquared + ...
                          CbDividedByCr);
    
    eyeMapLuma = im2double(zeros(size(grayImage)));
%     eyeMapLuma = im2double(grayImage);
%     grayImage = imadjust(grayImage);
    for n=1 : 10
        kernel = strel('disk', n);
%         kernel = strel('ball', n, n);
        dilation = imdilate(Y, kernel);
%         figure; imshow(dilation); title('dilation'); pause;
        erosion = imerode(Y, kernel);
%         figure; imshow(erosion); pause; title('erosion'); pause;
%         figure; imshow(im2double(dilation./(1 + erosion)));  pause;
        eyeMapLuma = eyeMapLuma + n*im2double(dilation./(1 + erosion));
%         eyeMapLuma = imadjust(eyeMapLuma);
    end
%     eyeMapLuma = imadjust(eyeMapLuma);
    
    grayFace = im2double(imadjust(rgb2gray(face)));
    filtered = im2double(zeros(size(grayFace)));
    for k = 1:10
        h = (ones(k, k) - fspecial('gaussian', k, 0.1))./k^2;
        filtered = filtered + im2double(imfilter(grayFace, h)) / k;
        filtered = imadjust(filtered);
%         kernel = strel('disk', k);
%         dilation = imdilate(grayFace, kernel);
%         erosion = imerode(grayFace, kernel);
%         filtered = filtered + im2double(erosion./(1 + dilation));
%         filtered = imadjust(filtered);
%         figure; imshow( filtered ); title('filtered'); pause;

    end
%     filtered = (filtered < 1) & faceMask;
    filtered = im2double(ones(size(filtered))) - filtered;
    filtered = filtered.*faceMask;
%     filtered = imdilate(, kernel);
%     figure; imshow( filtered ); title('filtered'); pause;
%     figure; imshow(eyeMapLuma); title('eyeMapLuma'); pause
    eyeMapLuma = imdilate(eyeMapLuma, strel('disk', 1));
%     eyeMapLuma = eyeMapLuma .* sqrt(filtered);
%     eyeMapLuma = eyeMapLuma .* filtered;
%      eyeMapLuma = imadjust(eyeMapLuma);
%      figure; imshow(eyeMapLuma); title('eyeMapLuma'); pause
     
    eyeMapChroma = histeq(eyeMapChroma);
%     figure; imshow(eyeMapChroma); title('eyeMapChroma'); pause;
%     eyeMapChroma = histeq(eyeMapChroma);
%     figure; imshow(eyeMapChroma); title('eyeMapChroma'); pause;
     
%     G = fspecial('gaussian',[5 5],2);
%     eyeMapLuma = imfilter(eyeMapLuma, G, 'same');

%     figure; imshow(eyeMapLuma); title('eyeMapLuma'); pause;
%     eyeMapLuma = histeq(eyeMapLuma);
%     eyeMapLuma = imadjust(eyeMapLuma);
%     figure; imshow(eyeMapLuma); title('eyeMapLuma'); pause;
    
    eyeMap = imadjust(im2double(eyeMapChroma) .* im2double(eyeMapLuma));
%     eyeMap = im2double(eyeMapChroma) .* im2double(eyeMapLuma);
    originalEyeMap = eyeMap;
%     figure; imshow(originalEyeMap); title('originalEyeMap'); pause;
    eyeMap = eyeMap.*filteredFaceMaskCopy;
    
    eyeMap = eyeMap.*im2double(faceMask);
%     figure; imshow(eyeMap); title('eyeMap'); pause;
    
%     eyeMap = histeq(im2double(eyeMapChroma) .* im2double(eyeMapLuma));
%     
%     eyeMap = eyeMap.*im2double(faceMask);
%     figure; imshow(eyeMap); title('eyeMap'); pause;


    CrSquared = Cr.^2;
    n = length(CrSquared(:));
    numerator = (1 / n) * sum(CrSquared(:));
    CrDividedByCb = Cr./Cb;
    denumerator = (1 / n) * sum(CrDividedByCb(:));
    nn = 0.95 * (numerator / denumerator)
    mouthMap = CrSquared.*(CrSquared - nn*CrDividedByCb).^2;
    mouthMap = imadjust(mouthMap);
%     figure; imshow(mouthMap); title('mouthMap'); pause;
%     mouthMap = or(imfilter(mouthMap, [-1 0 1]'), imfilter(mouthMap, [1 0 -1]'));
    
%     figure; imshow(mouthMap); title('mouthMap'); pause;
    mouthMap = mouthMap.*faceMask;
    mouthMap = mouthMap > 0.8 * max(mouthMap(:));
    mouthMap = imdilate(mouthMap, strel('disk', 4));
    mouthMap = imfill(mouthMap, 'holes');
    mouthMap = imerode(mouthMap, strel('disk', 4));
   
%     mouthMapRep = repmat(mouthMap, [1,1,3]);
%     mouthMap = rgb2gray(face) .* mouthMap;
%     mouthMap = or(imfilter(mouthMap, [-1 0 1]'), imfilter(mouthMap, [1 0 -1]'));
%     mouthMap = imdilate(mouthMap, strel('disk', 2));
%     mouthMap = imerode(mouthMap, strel('disk', 2));
%     figure; imshow(mouthMap); title('mouthMap'); pause;
    mouthMap = imfill(mouthMap, 'holes');
%     figure; imshow(mouthMap); title('mouthMap'); pause;
%     mouthMap = ExtractNLargestBlobs(mouthMap, 1);
    
%     mouthMap = imerode(mouthMap, strel('disk', 4));
%     figure; imshow(mouthMap); title('mouthMap'); pause;
    
%     figure; imshow(eyeMap); title('eyeMap'); pause;
    eyeMap = eyeMap - mouthMap;
%     eyeMap = eyeMap .* mouthMap;
%     figure; imshow(eyeMap); title('eyeMap - mouthMap'); pause;
    
    
%     figure; imshow(eyeMap); title('eyeMap11'); pause;
%     eyeMap = eyeMap - and(overSaturatedMask, eyeMap);
%     figure; imshow(eyeMap); title('eyeMap2'); pause;
%     eyeMap = imfill(eyeMap, 'holes');
%     eyeMap = imdilate(eyeMap, strel('disk', 4));
%     eyeMap = ExtractNLargestBlobs(eyeMap, 2);
%     figure; imshow(eyeMap); title('final Eyemap'); pause;
%     figure; imshow(eyeMap); title('final Eyemap'); pause;
%     eyeMap = imerode(eyeMap, strel('disk', 3));
%     figure; imshow(eyeMap); title('final Eyemap'); pause;
%     eyeMap = imdilate(eyeMap, strel('disk', 14));
%     figure; imshow(eyeMap); title('final Eyemap'); pause;
%     eyeMap = imfill(eyeMap, 'holes');
%     figure; imshow(eyeMap); title('final Eyemap'); pause;
%     eyeMap = ExtractNLargestBlobs(eyeMap, 2);
%     figure; imshow(eyeMap); title('final Eyemap'); pause;
%     
%     eyeMapRep = repmat(eyeMap, [1,1,3]);
%     I = face.*uint8(eyeMapRep);
    
%     figure; imshow(estimatedSkinMask); title('estimatedSkinMask'); pause;
%     estimatedSkinMask = imerode(estimatedSkinMask, strel('disk', 5));
%     figure; imshow(estimatedSkinMask); title('estimatedSkinMask'); pause;
%     figure; imshow(eyeMap); title('eyeMap'); pause;
    estimatedSkinMask = imerode(estimatedSkinMask, strel('disk', 5));
    eyeMap = eyeMap .* ~estimatedSkinMask;
%     figure; imshow(eyeMap); title('eyeMap'); pause;
%     figure; imshow(eyeMap); title('eyeMap'); pause;
    eyeMap = imfill(eyeMap, 'holes');
%     figure; imshow(eyeMap); title('eyeMap'); pause;
%     figure; imshow(eyeMap); title('eyeMap'); pause;
    eyeMap = eyeMap > 0;
%     figure; imshow(eyeMap); title('eyeMap'); pause;
%     figure; imshow(eyeMap); title('eyeMap'); pause;
    
    eyeMap = imdilate(eyeMap, strel('disk', 10));
%     figure; imshow(eyeMap); title('eyeMap'); pause;
% %     figure; imshow(eyeMap); title('eyeMap'); pause;
%     eyeMap = eyeMap >= 0.95 * max(eyeMap(:));
%     eyeMap = and(originalEyeMap, eyeMap);
    
%     eyeMapMaskZero = eyeMap == 0;
%     figure; imshow(faceMask); title('faceMask'); pause;
    

%     faceHSV = rgb2ycbcr(face);
%     faceH = faceHSV(:, :, 1);
%     meanFaceH = mean(faceH(faceMask));
%     
%     faceS = faceHSV(:, :, 2);
%     meanFaceS = mean(faceS(faceMask));
%    
%     skinH = (faceH < meanFaceH + 35) & ...
%             (faceH > meanFaceH - 35) & ... 
%             (faceS < meanFaceS + 10) & ...
%             (faceS > meanFaceS - 10) ;
%     figure; imshow( skinH ); title('skinH'); pause;
%     faceHSV = rgb2hsv(face);
%     R = face(:,:,1);
%     G = face(:,:,1);
%     B = face(:,:,1);
%     skinH = (R + G + B) < 170 & faceMask;
%     figure; imshow( skinH ); title('skinH'); pause;
%     skinH = imdilate(skinH , strel('disk', 10));
%     figure; imshow( skinH ); title('skinH'); pause;
    
%     skinH = imfill(skinH, 'holes');
    

    
    
    
%     cc = and(im2double(eyeMapChroma), im2double(eyeMapLuma));
% figure; imshow(eyeMap); title('eyeMap'); pause;
    eyeMapRep = repmat(eyeMap, [1,1,3]);
    I = rgbImage.*uint8(eyeMapRep);
%     I = rgbImage;
% figure; imshow(originalEyeMap); title('originalEyeMap'); pause;
    IGradient = imgradient(originalEyeMap);
%      figure; imshow(IGradient); title('IGradient'); pause;
    IGradient = IGradient.*eyeMap;
%      figure; imshow(IGradient); title('IGradient'); pause;
    IGradient = IGradient.*filteredFaceMaskEyes;
%     figure; imshow(IGradient); title('IGradient'); pause;
%     IGradient = IGradient.*skinH;
%     IGradient = IGradient.*filtered;
%     faceMask = imerode(faceMask, strel('disk',4));
    IGradient = IGradient.*faceMask;
%     nonSkinMask1 = Cb >= Cr;
%     nonSkinMask1 = imdilate(nonSkinMask1, strel('disk', 10));
%     IGradient = IGradient.*nonSkinMask1;
%     IGradient = IGradient.*meanColorizedFace;
%     figure; imshow(nonSkinMask1); title('nonSkinMask1'); pause;
%     figure; imshow(I); title('I dark'); pause;
    
    IBright = IGradient;
    IDark = IGradient;
%     
%     [centersDark, radiiDark] = imfindcircles(IBright, [4 40], ...
%                  'ObjectPolarity', 'dark', 'sensitivity', 0.99);
%    
%     
%     mask = originalEyeMap > 0.2 *(max(originalEyeMap(:)));
%     eyeMap = eyeMap .* mask;
%     eyeMap = imdilate(eyeMap, strel('disk', 10));
%     eyeMapRep = repmat(eyeMap, [1,1,3]);
%     IBright = rgbImage.*uint8(eyeMapRep);
    
%     originalEyeMap = originalEyeMap .* eyeMap;
%     IBright = originalEyeMap > 0.5 * max(originalEyeMap(:));
%     IBright = imdilate(IBright, strel('disk', 20));
%     IBright = imfill(IBright, 'holes');
%     eyeMap = IBright;
% %         eyeMap = imdilate(eyeMap, strel('disk', 10));
%     eyeMapRep = repmat(eyeMap, [1,1,3]);
%     IBright = rgbImage.*uint8(eyeMapRep);
     figure; imshow(IBright); title('IBright'); pause;

    [centersBright, radiiBright] = imfindcircles(IBright, [4 20], ...
                 'ObjectPolarity', 'bright', 'sensitivity', 0.99);
             
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
    
    figure; imshow(face); title('face before rotated'); 
    viscircles([leftEyeCenterX, leftEyeCenterY], leftRadius,'EdgeColor','r');
    viscircles([rightEyeCenterX, rightEyeCenterY], rightRadius,'EdgeColor', 'g');
    pause;
    
    referenceDirection = [0 1];
    eyeDirection = [leftEyeCenterX-rightEyeCenterX, leftEyeCenterY-rightEyeCenterY];
    eyeDirection = eyeDirection / norm(eyeDirection);
    
    angle = acos(dot(referenceDirection, eyeDirection)) * 180 / pi;
    angle = 90 - angle;

%     figure; imshow(face); title('face before rotated'); pause;
    face = imrotate(face, -angle);
%     figure; imshow(face); title('face after rotated'); pause;
    
    output = rgbImage;
end
