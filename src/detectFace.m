function output = detectFace(rgbImage)
%     figure; imshow(rgbImage); title('rgbImage');    
    
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
    
    % Orginal 0.47
    
%     faceMask2 = imdilate(faceMask, strel('disk',10));
%     faceMaskRep2 = repmat(faceMask2, [1,1,3]);
%     face2 = rgbImage.*uint8(faceMaskRep2);
%     

%     grayFaceEyes = imadjust(rgb2gray(face));
% %     figure; imshow(grayFace); title('grayFace'); pause;
%     grayFaceEyes = imfilter(grayFaceEyes, fspecial('average', 2));
% %     figure; imshow(grayFace); title('grayFace'); pause;
% %     grayFace = imfilter(grayFace, fspecial('gaussian', 5));
% %     figure; imshow(grayFace); title('grayFace'); pause;
% %     grayFace = imfilter(grayFace, fspecial('gaussian'));
% %     figure; imshow(grayFace); title('grayFace'); pause;
% %     grayFace = imfilter(grayFace, fspecial('gaussian'));
% %     figure; imshow(grayFace); title('grayFace'); pause;
%     grayFaceEyes = im2bw(grayFaceEyes, 0.47);
%     filtEyes = imfilter(grayFaceEyes, fspecial('laplacian'));
%     figure; imshow(filtEyes); title('filtEyes'); pause;
%     [labeledImage, numberOfBlobs] = bwlabel(filtEyes);
%     convexProperties = regionprops(labeledImage, 'ConvexArea');
%     convexAreas = [convexProperties.ConvexArea];
%     [sortedConvexAreas, sortedIndices] = sort(convexAreas, 2, 'descend');
%     filtEyesMask = ismember(labeledImage, sortedIndices(1));
%     filtEyes = filtEyes - filtEyesMask;
%    
%     filtEyes = imfill(filtEyes, 'holes');
% %     filtEyes = imdilate(filtEyes, strel('disk', 10));
% %     filtEyes = bwareaopen(filtEyes, 30);
% %     filtEyes = imdilate(filtEyes, strel('disk', 10));
% %     filtEyes = imfill(filtEyes, 'holes');
%     figure; imshow(filtEyes); title('filtEyes'); pause;
    
    grayFace = imadjust(rgb2gray(face));
%     figure; imshow(grayFace); title('grayFace'); pause;
%     grayFace = imfilter(grayFace, fspecial('average',2));
%     figure; imshow(grayFace); title('grayFace'); pause;
%     grayFace = imfilter(grayFace, fspecial('gaussian', 5));
%     figure; imshow(grayFace); title('grayFace'); pause;
%     grayFace = imfilter(grayFace, fspecial('gaussian'));
%     figure; imshow(grayFace); title('grayFace'); pause;
%     grayFace = imfilter(grayFace, fspecial('gaussian'));
%     figure; imshow(grayFace); title('grayFace'); pause;
    grayFace = im2bw(grayFace, 0.47);
    filteredFaceMask = imfilter(grayFace, fspecial('laplacian'));
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
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMaskFilled = imfill(filteredFaceMask, 'holes');
%     figure; imshow(filteredFaceMaskFilled); title('filteredFaceMaskFilled'); pause;
    filteredFaceMaskFilledEroded = imerode(filteredFaceMaskFilled, strel('disk', 1));
%     figure; imshow(filteredFaceMaskFilledEroded); title('filteredFaceMaskFilledEroded'); pause;
    filteredFaceMaskBorder = xor(filteredFaceMaskFilled, filteredFaceMaskFilledEroded);
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    
    originalFilteredFaceMask = imdilate(originalFilteredFaceMask, strel('disk', 1));
    originalFilteredFaceMask = imerode(originalFilteredFaceMask, strel('disk', 1));
    filteredFaceMaskCopy2 = originalFilteredFaceMask - and(originalFilteredFaceMask, filteredFaceMaskBorder); 
    filteredFaceMaskCopy2 = imfill(filteredFaceMaskCopy2, 'holes');
    
    filteredFaceMaskEyes = filteredFaceMask;
    filteredFaceMaskEyes = imfill(filteredFaceMaskEyes, 'holes');
    
    
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
    
%     skinMask1 = trueColorize(face, faceMask, [0.5, 0.5, 0.5]);
% %     skinMask1 = imdilate(skinMask1, strel('disk', 2));
% %     skinMask1 = imerode(skinMask1, strel('disk', 5));
%     figure; imshow(skinMask1, []); title('skinMask1'); pause;
%     
%     eyesMask1 = imfill(skinMask1, 'holes');
% %     figure; imshow(eyesMask1, []); title('eyesMask1'); pause;
%     eyesMask1 = xor(eyesMask1, skinMask1);
%     figure; imshow(eyesMask1, []); title('eyesMask1'); pause;
% %     eyesMask1 = imerode(eyesMask1, strel('disk', 2));
%     eyesMask1 = ExtractNLargestBlobs(eyesMask1, 2);
% %     eyesMask1 = imdilate(eyesMask1, strel('disk', 20));
% %     eyesMask1 = bwareaopen(eyesMask1, 30);
%     figure; imshow(eyesMask1, []); title('eyesMask1'); pause;
















    
    
    

    
    
    
    
    
  
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
    
%     grayFace = im2double(imadjust(rgb2gray(face)));
%     filtered = im2double(zeros(size(grayFace)));
%     for k = 1:10
%         h = (ones(k, k) - fspecial('gaussian', k, 0.1))./k^2;
%         filtered = filtered + im2double(imfilter(grayFace, h)) / k;
%         filtered = imadjust(filtered);
%     end
%     filtered = im2double(ones(size(filtered))) - filtered;
%     filtered = filtered.*faceMask;
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
    
%     eyeMap = eyeMap.*filtEyes;
%     figure; imshow(skinMask1, []); title('skinMask1'); pause;
%     figure; imshow(~skinMask1, []); title('~skinMask1'); pause;
%     figure; imshow(eyeMap, []); title('eyeMap'); pause;
%     eyeMap = eyeMap.*(eyesMask1);
%     figure; imshow(eyeMap, []); title('eyeMap.*~skinMask1'); pause;
    
    eyeMap = eyeMap.*filteredFaceMaskCopy;
%     figure; imshow(filteredFaceMaskCopy); title('filteredFaceMaskCopy'); pause;
%     figure; imshow(eyeMap); title('eyeMap.*filteredFaceMaskCopy'); pause;
    
    eyeMap = eyeMap.*filteredFaceMaskCopy2;
%     figure; imshow(filteredFaceMaskCopy2); title('filteredFaceMaskCopy2'); pause;
%     figure; imshow(eyeMap); title('eyeMap.*filteredFaceMaskCopy2'); pause;
    
    eyeMap = eyeMap.*im2double(faceMask);
%     figure; imshow(eyeMap); title('eyeMap.*im2double(faceMask)'); pause;
    



    n = length(nonzeros(faceMask));
    
    Cb = double(Cb_orig);
    Cr = double(Cr_orig);
    
    CrSquared = Cr.^2;
    numerator = (1 / n) * sum(CrSquared(:));
    
    CrDividedByCb = Cr./Cb;
    denumerator = (1 / n) * sum(CrDividedByCb(:));
    
    nn = 0.95 * (numerator / denumerator);
    
    mouthMap = CrSquared.*(CrSquared - nn*CrDividedByCb).^2;
    
    mouthMap = im2double(mouthMap);
    mouthMap = 255 * ((mouthMap - min(mouthMap(:))) ...
                 ./ (max(mouthMap(:)) - min(mouthMap(:))));
    mouthMap = uint8(mouthMap);
    originalMouthMap = mouthMap;
    
    
    
    
    
    
%     
%     mouthMask = mouthMap > 0.5 * max(mouthMap(faceMask));
% %     figure; imshow(mouthMask); title('mouthMask'); pause;
%     [labeledImage, numberOfBlobs] = bwlabel(faceMask);
%     bb = regionprops(labeledImage, 'BoundingBox');
% %     bb(1).BoundingBox(1)
% %     bb(1).BoundingBox(2)
%     width = bb(1).BoundingBox(2);
%     mouthMask = imdilate(mouthMask, strel('disk', round(width/4)));
% %     figure; imshow(mouthMask); title('mouthMask'); pause;
%     figure; imshow(eyeMap); title('eyeMap'); pause;
% %     eyeMap = eyeMap.*~mouthMask;
%     figure; imshow(eyeMap); title('eyeMap'); pause;
    
    
%     figure; imshow(eyeMap); title('eyeMap11'); pause;
%     eyeMap = eyeMap - and(overSaturatedMask, eyeMap);
    eyeMap = eyeMap.*~overSaturatedMask;
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
    
%     figure; imshow(face); title('face');
%     viscircles([leftEyeCenterX, leftEyeCenterY], leftRadius,'EdgeColor','r');
%     viscircles([rightEyeCenterX, rightEyeCenterY], rightRadius,'EdgeColor', 'g');
% %     viscircles([mouthX, mouthY], rightRadius,'EdgeColor', 'b');
%     pause;
    
    eyeDistance = rightEyeCenterX - leftEyeCenterX;
    measureDistance = eyeDistance * 0.8;
    eyeY = (leftEyeCenterY + rightEyeCenterY) / 2;

    eyeRegionMask = zeros(size(face(:,:,1)));
    eyeRegionMask(1: eyeY+measureDistance, :) = 1;

    mouthMap = mouthMap .* uint8(faceMask);

    mouthMap = mouthMap > 0.5 * max(mouthMap(~eyeRegionMask));

    mouthMap = imdilate(mouthMap, strel('disk', 4));
    mouthMap = imfill(mouthMap, 'holes');
    mouthMap = imerode(mouthMap, strel('disk', 4));



    

    BW = and(mouthMap, ~eyeRegionMask);
    
    BW = ExtractNLargestBlobs(BW, 1);
    s = regionprops(BW, mouthMap, {'Centroid','WeightedCentroid', 'Area'});

    mouthX = s(1).Centroid(1);
    mouthY = s(1).Centroid(2);

    
    
    
    
    
    
%     figure; imshow(face); title('face before rotated'); 
%     viscircles([leftEyeCenterX, leftEyeCenterY], leftRadius,'EdgeColor','r');
%     viscircles([rightEyeCenterX, rightEyeCenterY], rightRadius,'EdgeColor', 'g');
%     viscircles([mouthX, mouthY], rightRadius,'EdgeColor', 'b');
%     pause;
    
    referenceDirection = [0 1];
    eyeDirection = [leftEyeCenterX-rightEyeCenterX, leftEyeCenterY-rightEyeCenterY];
    eyeDirection = eyeDirection / norm(eyeDirection);
    
    angle = acos(dot(referenceDirection, eyeDirection)) * 180 / pi;
    angle = 90 - angle;

%     figure; imshow(face); title('face before rotated'); pause;
    
    
    marker=zeros(size(face(:,:,1)));
    originalMarker = marker;
    
    face = imrotate(face, -angle);
    rgbImage = imrotate(rgbImage, -angle);
    
    marker(round(leftEyeCenterY), round(leftEyeCenterX))=1;
    marker_rot = imrotate(marker, -angle);
    [leftEyeCenterY, leftEyeCenterX]=find(marker_rot);
  
    marker = originalMarker;
    marker(round(rightEyeCenterY), round(rightEyeCenterX))=1;
    marker_rot = imrotate(marker, -angle);
    [rightEyeCenterY, rightEyeCenterX]=find(marker_rot);
    
    marker = originalMarker;
    marker(round(mouthY), round(mouthX))=1;
    marker_rot = imrotate(marker, -angle);
    [mouthY, mouthX]=find(marker_rot);

%     figure; imshow(face); title('face after rotated'); pause;
%     viscircles([leftEyeCenterX, leftEyeCenterY], leftRadius,'EdgeColor','r');
%     viscircles([rightEyeCenterX, rightEyeCenterY], rightRadius,'EdgeColor', 'g');
%     viscircles([mouthX, mouthY], rightRadius,'EdgeColor', 'b');
%     pause;
    
    
    eyeDistance = rightEyeCenterX - leftEyeCenterX;
    eyeY = (leftEyeCenterY + rightEyeCenterY) / 2;
    
    eyeMouthDistance = mouthY - eyeY;
    
    offsetX = eyeDistance * 0.3;
    offsetY = eyeMouthDistance * 0.2;
    
    output = rgbImage(round(eyeY-offsetY) : round(mouthY+offsetY), ...
                      round(leftEyeCenterX-offsetX) :  round(rightEyeCenterX+offsetX), :);
%     figure; imshow(output); title('output');
%     pause;
    
end
