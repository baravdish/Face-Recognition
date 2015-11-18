function output = detectFace(rgbImage)
%     figure; imshow(rgbImage); title('rgbImage'); pause;   
    
    rgbImage = padarray(rgbImage,[5 5],'both');
    
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
    figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    
    filteredFaceMask = bwareaopen(filteredFaceMask, 50);
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
   
    [labeledImage, numberOfBlobs] = bwlabel(filteredFaceMask);
    convexProperties = regionprops(labeledImage, 'ConvexArea');
    convexAreas = [convexProperties.ConvexArea];
    [sortedConvexAreas, sortedIndices] = sort(convexAreas, 2, 'descend');
    
    filteredFaceMask = ismember(labeledImage, sortedIndices(1));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMaskCopy = filteredFaceMaskCopy - filteredFaceMask;
    filteredFaceMaskCopy = imfill(filteredFaceMaskCopy, 'holes');
    filteredFaceMaskCopy = imdilate(filteredFaceMaskCopy, strel('disk', 5));
    filteredFaceMaskCopy = imfill(filteredFaceMaskCopy, 'holes');
    
    figure; imshow(filteredFaceMaskCopy); title('filteredFaceMaskCopy'); pause;
    
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk', 30));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask4'); pause;

    filteredFaceMask2 = filteredFaceMask;
    filteredFaceMask2 = imfill(filteredFaceMask2, 'holes');
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

    figure; imshow(face); title('face'); pause;
  
    Y = im2double(Y);
    Cr = im2double(Cr);
    Cb = im2double(Cb);
    
    CbSquared = Cb.^2;
    CrInvertedSquared = (1-Cr).^2;
    CbSquaredDiviedByCr = CbSquared./Cr;
    eyeMapChroma = (1/3)*(CbSquared + ...
                          CrInvertedSquared + ...
                          CbSquaredDiviedByCr);

    eyeMapLuma = im2double(zeros(size(grayImage)));
%     eyeMapLuma = im2double(grayImage);
    grayImage = imadjust(grayImage);
    for n=1 : 10
        kernel = strel('disk', n);
        dilation = imdilate(Y, kernel);
%         figure; imshow(dilation); title('dilation'); pause;
        erosion = imerode(Y, kernel);
%         figure; imshow(erosion); title('erosion'); pause;
        eyeMapLuma = eyeMapLuma + im2double(dilation./(1 + erosion));
    end
    
   
%     figure; imshow(eyeMapChroma); title('eyeMapChroma'); pause;
    eyeMapChroma = histeq(eyeMapChroma);
    figure; imshow(eyeMapChroma); title('eyeMapChroma'); pause;
     
%     G = fspecial('gaussian',[5 5],2);
%     eyeMapLuma = imfilter(eyeMapLuma, G, 'same');
    figure; imshow(eyeMapLuma); title('eyeMapLuma'); pause;
%     eyeMapLuma = histeq(eyeMapLuma);
%     ftigure; imshow(eyeMapLuma); title('eyeMapLuma'); pause;
    
    eyeMap = imadjust(im2double(eyeMapChroma) .* im2double(eyeMapLuma));
    eyeMap = eyeMap.*filteredFaceMaskCopy;
    
    eyeMap = eyeMap.*im2double(faceMask);
    figure; imshow(eyeMap); title('eyeMap'); pause;
    
%     eyeMap = histeq(im2double(eyeMapChroma) .* im2double(eyeMapLuma));
%     
%     eyeMap = eyeMap.*im2double(faceMask);
%     figure; imshow(eyeMap); title('eyeMap'); pause;



    n = length(CbSquared(:));
    numerator = (1 / n) * sum(CbSquared(:));
    CrDividedByCb = Cr./Cb;
    denumerator = (1 / n) * sum(CrDividedByCb(:));
    nn = 0.95 * (numerator / denumerator)
    mouthMap = CbSquared.*(CbSquared - nn*CrDividedByCb).^2;
    mouthMap = imadjust(mouthMap);
    figure; imshow(mouthMap); title('mouthMap'); pause;

    
%     eyeMapRep = repmat(eyeMap, [1,1,3]);
%     I = face.*uint8(eyeMapRep);
    
    I = eyeMap;
    [centersDark, radiiDark, metric] = imfindcircles(I, [4 40], ...
                 'ObjectPolarity', 'bright', 'sensitivity', 0.99);

    figure; imshow(I); title('I'); pause;
%     viscircles(centersDark(1:2,:), radiiDark(1:2),'LineStyle','--');
%     viscircles(centersDark, radiiDark,'LineStyle','--');
%     pause;
    
     
    eyeMapSize = size(I)
    imageSizeX = eyeMapSize(1,1,1)
    imageSizeY = eyeMapSize(1,2,1)
    [cols rows ] = meshgrid(1:imageSizeY, 1:imageSizeX);

    centerX1 = centersDark(1,1,1)
    centerY1 = centersDark(1,2,1)
    radius1 = radiiDark(1,1)
    irisMask = ((rows - centerY1).^2 + (cols - centerX1).^2) <= radius1.^2;
    
    centerX2 = centersDark(2,1,1)
    centerY2 = centersDark(2,2,1)
    radius2 = radiiDark(2,1)
    
    for n=2 : size(centersDark, 1)
        centerX2 = centersDark(n,1,1)
        centerY2 = centersDark(n,2,1)  
        radius2 = radiiDark(n,1) 
        
        distance = sqrt((centerX2-centerX1)^2 + (centerY2-centerY1)^2);
        
        if distance > (radius1 + radius2)*2
            n = size(centersDark, 1)
            break
        end
    end
         

    irisMask = irisMask | ((rows - centerY2).^2 + (cols - centerX2).^2) <= radius2.^2;
    
    irisMaskRep = repmat(irisMask, [1,1,3]);
    irisMap = face.*uint8(irisMaskRep);
    
    figure; imshow(irisMap); title('irisMap'); pause;
    
    
    
    output = rgbImage;
end
