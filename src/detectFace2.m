    foregroundMask = detractBackground(input);
    
    faceMask = and(faceMask, foregroundMask);
       
    imshow(faceMask); title('faceMask'); pause;

    lapImg = imfilter(imadjust(rgb2gray(input)), fspecial('laplacian'));

    imshow(lapImg); title('lapImg'); pause;

    figure
    level = graythresh(input);
    bw = im2bw(input, level);
    bw = bwareaopen(bw, 50);
    bw = ~bw;
    bw = imfill(bw, 'holes');
    imshow(bw)
    pause
    
    foreground = bw;
    background = ~bw;
    
    nPixels = size(lapImg, 1) * size(lapImg, 2);

    rega = lapImg.*uint8(bw);
    inva = lapImg.*uint8(~bw);

    rega = sum(rega(:)) / length(nonzeros(rega(:)))
    inva = sum(inva(:)) / length(nonzeros(inva(:)))
    
    homogeneousBackground = 0;
    
    lowa = rega;
    homo = inva
    
    if inva > rega
       foreground = ~bw; 
       background = bw;
       homo = rega
       lowa = inva;
    end
    
    if homo < 8 && lowa / homo  > 3
        homogeneousBackground = 1;
    end
    
    foreground = imdilate(foreground, strel('disk', 10));
    foreground = imfill(foreground, 'holes');
    foreground = imerode(foreground, strel('disk', 10));
    
    figure; imshow(foreground); title('foreground'); pause;
    
    figure; imshow(faceMask); title('faceMask'); pause;
    faceMask = imfill(faceMask, 'holes');
    
    faceMask = ExtractNLargestBlobs(faceMask, 1);

    figure; imshow(faceMask); title('faceMask'); pause;

    faceMask = imdilate(faceMask, strel('disk', 20));
    faceMask = imerode(faceMask, strel('disk', 25));
    
    if homogeneousBackground == 1
        faceMask = and(foreground, faceMask);
    end
    
    imshow(faceMask); title('faceMask'); pause;

    faceMaskCopy = faceMask;
    minCol = [];
    maxCol = [];
    minRow = [];
    maxRow = [];
    val = 80;
    while( isempty(minCol) || isempty(maxCol) || ... 
           isempty(minRow) || isempty(maxRow) )
        faceMaskCopy = faceMask;
        minCol = [];
        maxCol = [];
        minRow = [];
        maxRow = [];
        faceMaskCopy = imerode(faceMaskCopy, strel('disk', val));
        faceMaskCopy = imdilate(faceMaskCopy, strel('disk', val));
        [row, col] = find(faceMaskCopy(:,:,1) ~= 0);
        minCol = min(col);
        maxCol = max(col);
        minRow = min(row);
        maxRow = max(row);
        val = max(val - 10, 1);
%         imshow(faceMaskCopy); title('faceMaskCopy'); pause;
    end
    
    faceMask = faceMaskCopy;

    faceMask = ExtractNLargestBlobs(faceMask, 1);
    
    figure; imshow(faceMask); title('cropped faceMask'); pause;
    figure; imshow(faceMask); title('faceMask');  pause;
    
    newMask = meanColorize(input, faceMask);

    newMask = and(newMask, faceMask);
    
    figure; imshow(newMask); title('newMask'); pause;
    newMask = ExtractNLargestBlobs(newMask, 1);
    newMask = imdilate(newMask, strel('disk', 4));
    newMask = imfill(newMask, 'holes');
    
    newMask = imerode(newMask, strel('disk', 7));
    figure; imshow(newMask); title('newMask morphed'); pause;
    faceMask = and(faceMask, newMask);
    figure; imshow(faceMask); title('faceMask');  pause;
    finalFaceMask = faceMask;
    faceMask = repmat(faceMask, [1,1,3]);
    face = input.*uint8(faceMask);
    
    figure; imshow(face); title('face');
    
    filteredFaceMask = imfilter(im2bw( imadjust(rgb2gray(face))), fspecial('laplacian'));
    filteredFaceMask = ExtractNLargestBlobs(filteredFaceMask, 1);
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',30));
    filteredFaceMask = ~filteredFaceMask;
    filteredFaceMask = and(filteredFaceMask, im2bw(face) > 0);
    filteredFaceMask = ExtractNLargestBlobs(filteredFaceMask, 1);
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));

    mask = filteredFaceMask;

    faceMask = and(mask, faceMask(:,:,1));
    redChannel = face(:, :, 1);
    greenChannel = face(:, :, 2);
    blueChannel = face(:, :, 3);
   
    faceMask = ExtractNLargestBlobs(faceMask, 1);
    void = ExtractNLargestBlobs(~faceMask, 1);
    faceMask = or(~void, faceMask);
    maskedRed = uint8(redChannel) .*uint8(faceMask);
    maskedGreen = uint8(greenChannel) .*uint8(faceMask);
    maskedBlue = uint8(blueChannel) .*uint8(faceMask);
    maskedRgbImage = cat(3, maskedRed, maskedGreen, maskedBlue);
    face = maskedRgbImage;
    
    faceHSV = rgb2hsv(maskedRgbImage);
    H = faceHSV(:,:,1);
    S = faceHSV(:,:,2);
    V = faceHSV(:,:,3);

    dominantHsv = [round(255*mean2(faceHSV(:,:,1))) ...
                   round(255*mean2(faceHSV(:,:,2))) ...
                   round(255*mean2(faceHSV(:,:,3)))];
    dominantHsv = (dominantHsv - min(dominantHsv))./(max(dominantHsv) - min(dominantHsv));
    figure
    imshow(faceHSV);
    title('faceHSV');

    [row, col] = find(face(:,:,1) ~= 0);
   
    minCol = min(col);
    maxCol = max(col);
    
    minRow = min(row);
    maxRow = max(row);
    
    cropImg = input(minRow:maxRow, minCol:maxCol, :);
    faceCropped = face(minRow:maxRow, minCol:maxCol, :);
    
    figure;
    imshow(face);
    title('face');
    pause;

    
    i=input;
    iMask = face > 0;
    iMask = or(or( iMask(:,:,1), iMask(:,:,2)), iMask(:,:,3));
    iMaskRep = repmat(iMask, [1,1,3]);

    figure;
    subplot(4,4,1)
    imshow(i)
    title('original image');
    iycbcr=rgb2ycbcr(i);
    iycbcr = im2double(iycbcr);
    subplot(4,4,2)
    imshow(iycbcr)
    title('YCBCR space');
    y=iycbcr(:,:,1);
    cb=iycbcr(:,:,2);
    cr=iycbcr(:,:,3);
    ccb=cb.^2;
    subplot(4,4,3)
    imshow(ccb)
    title('CB^2');
    ccr=(1-cr).^2;
    subplot(4,4,4)
    imshow(ccr)
    title('(1-CR)^2');
    cbcr=ccb./cr;
    subplot(4,4,5)
    imshow(cbcr)
    title('CB/CR');
    igray=rgb2gray(i);
    subplot(4,4,6)
    imshow(igray)
    title('Gray space');
    g=1./3;
    l=g*ccb;
    m=g*ccr;
    n=g*cbcr;

    eyemapchr=l+m+n;
    subplot(4,4,7)

    imshow(eyemapchr)
    title('Chrom Eyemap');
    J=histeq(eyemapchr);
    subplot(4,4,8)
    eyemapchr = J;

    imshow(J)
    title('Equalized image');
    
    SE=strel('disk',15,8);
    o=imdilate(igray,SE);
    p=1+imerode(igray,SE);
    eyemaplum=o./p;
    

    subplot(4,4,9)
    imshow(eyemaplum)
    title('Lum Eyemap');
    pause;
    
    figure
    imshow(eyemapchr); title('eyemapchr Eyemap'); pause;
    imshow(eyemaplum); title('eyemaplum Eyemap'); pause;

    
    vea = eyemapchr.*finalFaceMask;
    imshow(vea); title('vea'); pause;
    vec = nonzeros(vea(:));
    soretedVec = sort(vec, 'descend');
    eyemapchrTh = 0.95 * soretedVec(1,1)
    
    eyemapchr22 = imgradient(rgb2gray(input)) / 255;
    imshow(eyemapchr22);
    title('imgradient');
    pause;


    figure; imshow(eyemapchr); title('eyemapchr before'); pause;

    figure; imshow(eyemapchr); title('eyemapchr after'); pause;
    eyemapchr = eyemapchr > eyemapchrTh;
    figure; imshow(eyemapchr); title('eyemapchr'); pause;
    
    vea = eyemaplum.*uint8(finalFaceMask);
    imshow(vea); title('vea'); pause;
    vec = nonzeros(vea(:));
    soretedVec = sort(vec, 'descend');
    eyemaplumTh = 0.9 * (soretedVec(1,1) / 255)

    eyemaplum = eyemaplum > eyemaplumTh;
    figure; imshow(eyemaplum); title('eyemaplum'); pause;

    
    eyemapchr = and(eyemapchr, faceMask);
  
    
    eyemaplum = and(eyemaplum, faceMask);
    
    

    imshow(eyemapchr); title('eyemapchr true'); pause;
    imshow(eyemaplum); title('eyemaplum true'); pause;
    eyeMap = and(eyemapchr, eyemaplum);
    imshow(eyeMap); title('eyeMap'); pause;

    figure;
    overSaturated = rgb2gray(face);
    
    imshow(imfilter(imadjust(overSaturated), fspecial('laplacian'))); title('filter imadjust(overSaturated)'); pause;
    overSaturated = imfilter(overSaturated, fspecial('laplacian'));

    overSaturated = imdilate(overSaturated, strel('disk', 2));
    imshow(overSaturated); title('dilate Over saturation'); pause;
    overSaturated = ~overSaturated;
    imshow(overSaturated); title('filter Over saturation'); pause;
    
    overSaturated = and(finalFaceMask, overSaturated);
    imshow(overSaturated); title('Over saturation'); pause;
    
    imshow(eyeMap); title('eyeMap11'); pause;
    eyeMap = eyeMap - and(overSaturated,eyeMap);
    
    imshow(eyeMap); title('eyeMap2'); pause;

    figure;
    
    eyeMap = imfill(eyeMap, 'holes');

    eyeMap = imdilate(eyeMap, strel('disk', 4));

    eyeMap = ExtractNLargestBlobs(eyeMap, 2);

    imshow(eyeMap); title('final Eyemap'); pause;
    erodeKernel = strel('disk', 11);
    imshow(eyeMap); title('final Eyemap'); pause;
    eyeMap = imerode(eyeMap, strel('disk', 3));
    imshow(eyeMap); title('final Eyemap'); pause;
    eyeMap = imdilate(eyeMap, strel('disk', 14));
    imshow(eyeMap); title('final Eyemap'); pause;
    eyeMap = imfill(eyeMap, 'holes');
    imshow(eyeMap); title('final Eyemap'); pause;

    eyeMap = ExtractNLargestBlobs(eyeMap, 2);

    finalEyeMask = eyeMap;
    eyeMap = repmat(eyeMap, [1,1,3]);
    eyeMap = face.*uint8(eyeMap);
    imshow(eyeMap); title('final Eyemap'); 


    level = graythresh(eyeMap);
    I = im2bw(eyeMap, level);
    imshow(I); title('I Eyemap'); pause;

    Rmin=4; 
    Rmax=40;  
    [centersDark, radiiDark, metric] = imfindcircles(I, [Rmin Rmax], ...
                 'ObjectPolarity','dark','sensitivity',0.99);
     
    figure; imshow(eyeMap); title('final Eyemap'); pause;;
    viscircles(centersDark, radiiDark,'LineStyle','--');
    pause;
    
    eyeMapSize = size(eyeMap)
    imageSizeX = eyeMapSize(1,1,1)
    imageSizeY = eyeMapSize(1,2,1)
    [cols rows ] = meshgrid(1:imageSizeY, 1:imageSizeX);
%     irisMask = final;
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
    size(irisMask)
    
    irisMap = repmat(irisMask, [1,1,3]);
    irisMap = face.*uint8(irisMap);
    
    figure; imshow(irisMap); title('irisMap');
    
    eyeColor = [mean2(nonzeros(irisMap(:,:,1))) ...
                mean2(nonzeros(irisMap(:,:,2))) ...
                mean2(nonzeros(irisMap(:,:,3)))]
   
    eyeColor = rgb2hsv(eyeColor/255);

    referenceDirection = [1 0];
    eyeDirection = [centerX1-centerX2, centerY1-centerY2];
    eyeDirection = eyeDirection / norm(eyeDirection);
    angle = acos(dot(referenceDirection, eyeDirection)) * 180 / pi
    angle = min(angle, 180-angle);
    
    figure; imshow(face); title('face before rotated'); pause;
    face = imrotate(face, -angle);
    eyeMap = imrotate(eyeMap, -angle);
    irisMap = imrotate(irisMap, -angle);
    cropImg = imrotate(cropImg, -angle);
    faceCropped = imrotate(faceCropped, -angle);
    finalFaceMask = imrotate(finalFaceMask, -angle);
    finalEyeMask = imrotate(finalEyeMask, -angle);
    overSaturated = imrotate(overSaturated, -angle);
    figure; imshow(face); title('face after rotated'); pause;
    
    input = imrotate(input, -angle);
    iycbcr=rgb2ycbcr(input);
    iycbcr = im2double(iycbcr);
    y=iycbcr(:,:,1);
    cb=iycbcr(:,:,2);
    cr=iycbcr(:,:,3);
    

    cr2 = (cr).^2;
    n = length(cr2(:));
    ovan = (1 / n) * sum(cr2(:));
    cddcb = cr./cb;
    nedan = (1 / n) * sum(cddcb(:));
    
    figure
    nn = 0.95 * (ovan / nedan)
    mouthMap = cr2.*(cr2 - nn*cddcb).^2;
    
    mouthMap = imadjust(mouthMap);

    vea = mouthMap;
    imshow(vea); title('vea'); pause;
    vec = nonzeros(vea(:));
    soretedVec = sort(vec, 'descend');

    mouthTh = 0.99 * soretedVec(1,1)
    pause;
    imshow(mouthMap); title('mouthMap'); pause;
    mouthMap = mouthMap > mouthTh;

    imshow(mouthMap); title('mouthMap'); pause;

    mouthMap = or(imfilter(mouthMap, [-1 0 1]'), imfilter(mouthMap, [1 0 -1]'));
        
    mouthMap = uint8(mouthMap).*uint8(finalFaceMask);
    mouthMap = uint8(mouthMap) - uint8(mouthMap).*uint8(finalEyeMask);
    mouthMap = uint8(mouthMap) - uint8(mouthMap).*uint8(overSaturated);
    imshow(mouthMap); title('mouthMap'); pause;
    mouthMap = ExtractNLargestBlobs(mouthMap, 1);
    imshow(mouthMap); title('mouthMap'); pause;
    mouthMap = imdilate(mouthMap, strel('disk', 8));
    mouthMap = imfill(mouthMap, 'holes');
    mouthMap = repmat(mouthMap, [1,1,3]);
    mouthMap = face.*uint8(mouthMap);
    
    figure; imshow(mouthMap); title('mouthMap');
    
    eyesMouth = or(eyeMap(:,:,1), mouthMap(:,:,1));
    eyesMouth = repmat(eyesMouth, [1,1,3]);
    eyesMouth = face.*uint8(eyesMouth);
    
    figure; imshow(eyesMouth); title('eyesMouth');

    output = input;