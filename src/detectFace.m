%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: Image in which to detect face.
%
% output: Image containing only the face.
%%%%%%%%%%%%%%%%%%%%%%%%%%

% function output = deaatectFaceOrig(input)
%     
% %     figure; imshow(input); title('input'); 
%     
%     padding = 16;
%     paddedInput = padarray(input, [padding padding], 0, 'both');
%         
%     paddedInputYCbCr = rgb2ycbcr(paddedInput);
% %     paddedInputHsv = rgb2hsv(paddedInput);
% %     paddedInputLab = rgb2lab(paddedInput);
%     
% %     Y = paddedInputYCbCr(:,:,1); 
%     CB = paddedInputYCbCr(:,:,2); CR = paddedInputYCbCr(:,:,3);
%     
%     faceMaskWide = and(and(CB > 100, CB < 145), ... 
%                        and(CR > 136, CR < 165));
%     
%     faceMaskWideRep = repmat(faceMaskWide, [1,1,3]);
%     faceMaskWideImg = uint8(paddedInputYCbCr).*uint8(faceMaskWideRep);
% %     faceMaskWideImg(:,:,1) = imadjust(faceMaskWideImg(:,:,1));
% %     faceMaskWideImg(:,:,2) = imadjust(faceMaskWideImg(:,:,2));
% %     faceMaskWideImg(:,:,3) = imadjust(faceMaskWideImg(:,:,3));
%     
%     figure; imshow(faceMaskWideImg); title('faceMaskWideImg');
%     
%     meanColor = [mean2(nonzeros(faceMaskWideImg(:,:,1))) ...
%                  mean2(nonzeros(faceMaskWideImg(:,:,2))) ...
%                  mean2(nonzeros(faceMaskWideImg(:,:,3)))]
%  
%     redChannel = faceMaskWideImg(:,:,1);
%     greenChannel = faceMaskWideImg(:,:,2);
%     blueChannel = faceMaskWideImg(:,:,3);
%     
%     differnece = 8;
%     meanMask = (abs(redChannel - meanColor(1,1)) < differnece);
% %     meanMask = or(meanMask, abs(greenChannel - meanColor(1,2)) > differnece);
% %     meanMask = or(meanMask, (abs(blueChannel - meanColor(1,3)) > differnece));
%     skinMask = ~meanMask;
%     skinMask = imfill(skinMask, 'holes');
%     skinMaskOrig = skinMask;
%     figure; imshow(skinMask); title('skinMask'); 
%     
%     skinMask = imerode(skinMask, strel('disk', 3));
%     skinMask = imdilate(skinMask, strel('disk', 10));
%     skinMask = ExtractNLargestBlobs(skinMask, 1);
%     skinMask = imdilate(skinMask, strel('disk', 30));
%     skinMask = imfill(skinMask, 'holes');
%     
%      skinMask = imerode(skinMask, strel('disk', 30));
%      
%           skinMask = imerode(skinMask, strel('disk', 50));
%             skinMask = imdilate(skinMask, strel('disk', 50));
% %     skinMask = imerode(skinMask, strel('disk', 2));
% %     
% %     skinMask = imerode(skinMask, strel('disk', 30));
% %     skinMask = imdilate(skinMask, strel('disk', 30));
% %     skinMask = imdilate(skinMask, strel('disk', 7));
% %     skinMask = xor(skinMask, skinMaskOrig);
%     
%     figure; imshow(skinMask); title('skinMask updated'); 
%     
%     skinMaskRep = repmat(skinMask, [1,1,3]);
%     skinMaskImg = paddedInput.*uint8(skinMaskRep);
%     
%     figure; imshow(skinMaskImg); title('skinMaskImg'); 
%     
%     output = input;
% end

function output = detectFace(input)
<<<<<<< HEAD
    
%     figure
%     title('input');
%     imshow(input);
%     pause;
%     size(input)
    
%     input = create_padded_image(input, 0);
    input = padarray(input,[16 16],0,'both');
    
%     size(input)
%     title('input');
%     imshow(input);
%     pause;
%     
=======
  
    

>>>>>>> ae6fa4ab063fdb623d3b41a0b77a8379859ceaaf
    imgYCbCr = rgb2ycbcr(input);
    
    Y = imgYCbCr(:,:,1);
    CB = imgYCbCr(:,:,2);
    CR = imgYCbCr(:,:,3);
% 
%     faceMask = and(and(CB > 105, CB < 145), ... 
%                    and(CR > 136, CR < 165));
%                
    faceMask = and(and(CB > 100, CB < 145), ... 
                   and(CR > 136, CR < 165));
%     faceMask = and(and(CB > 100, CB < 145), ... 
%                    and(CR > 125, CR < 165));
        figure
    title('faceMask');
    imshow(faceMask);
%     pause;
    faceMask = imfill(faceMask, 'holes');
    
%     title('faceMask');
%     imshow(faceMask);
%     pause;
    faceMask = ExtractNLargestBlobs(faceMask, 1);
%     
%     faceMaskRep = repmat(faceMask, [1,1,3]);
%     input2 = input.*uint8(faceMaskRep);
    
%     input = lightNorm(input);
%       
%     title('input2');
%     imshow(input2);
%     pause;

%     erodeKernel = strel('disk', 30);
%     dilateKernel = strel('disk', 20);  
%     faceMask = imerode(faceMask, erodeKernel);
%     faceMask = imdilate(faceMask, dilateKernel);
    
%     for n = 1:6
%         erodeKernel = strel('disk', n);
%         dilateKernel = strel('disk', 1+n);
%         faceMask = imerode(faceMask, erodeKernel);
%         faceMask = imdilate(faceMask, dilateKernel);
%         
%         title('faceMask');
%         imshow(faceMask);
%         pause;
%     end
    
    figure
    title('faceMask');
    imshow(faceMask);
%     pause;
% 

    erodeKernel = strel('disk', 25);
    dilateKernel = strel('disk', 20);    
    faceMask = imdilate(faceMask, dilateKernel);
    faceMask = imerode(faceMask, erodeKernel);
    
%     title('faceMask');
%     imshow(faceMask);
%     pause;
    
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
        erodeKernel = strel('disk', val);
        dilateKernel = strel('disk', val); 
        faceMaskCopy = imerode(faceMaskCopy, erodeKernel);
        faceMaskCopy = imdilate(faceMaskCopy, dilateKernel);
        [row, col] = find(faceMaskCopy(:,:,1) ~= 0);
        minCol = min(col);
        maxCol = max(col);
        minRow = min(row);
        maxRow = max(row);
        val = max(val - 10, 1);

%         title('faceMaskCopy');
%         imshow(faceMaskCopy);
%         pause;
    end
    
    faceMask = faceMaskCopy;

    faceMask = ExtractNLargestBlobs(faceMask, 1);
    
%     figure
%     title('faceMask');
%     imshow(faceMask);
%     pause;
%     
%     figure;
%     title('faceMask');
%     imshow(faceMask);
    
    figure
    imshow(faceMask)
    faceMask = repmat(faceMask, [1,1,3]);
    face = input.*uint8(faceMask);
<<<<<<< HEAD
    
=======
        imshow(face)
    % Crop face
    [row, col] = find(face(:,:,1) ~= 0);
>>>>>>> ae6fa4ab063fdb623d3b41a0b77a8379859ceaaf
    
        figure;
    title('face');
    imshow(face);

    figure
%     imshow(face); title('face'); pause;
%     imshow(im2bw(face)); title('face'); pause;
    filteredFaceMask = imfilter(im2bw( imadjust(rgb2gray(face))), fspecial('laplacian'));
% filteredFaceMask = imfilter(filteredFaceMask, fspecial('laplacian'));
%     filteredFaceMask = imfilter(filteredFaceMask, fspecial('laplacian'));
    %     filteredFaceMask = imerode(filteredFaceMask, strel('disk',1));
%     imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMask = ExtractNLargestBlobs(filteredFaceMask, 1);
%     imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',30));
%     imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMask = ~filteredFaceMask;
%     imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMask = and(filteredFaceMask, im2bw(face) > 0);
%     imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMask = ExtractNLargestBlobs(filteredFaceMask, 1);
%     imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));

%      maska2 = imfill(im2bw(face), 'holes');
%      maska2 = imdilate(maska2, strel('disk',5));
%      maska2 = imfill(maska2, 'holes');
%      filteredFaceMask = and(filteredFaceMask, maska2);
%      filteredFaceMask = imfill(filteredFaceMask, 'holes');

%     filteredFaceMask = ExtractNLargestBlobs(filteredFaceMask, 1);
%     filteredFaceMask = imerode(filteredFaceMask, strel('disk',2));
%         filteredFaceMask = imfill(filteredFaceMask, 'holes');

% lbl = bwlabel(~BW);
% holes = ~(BW|ismember(lbl,unique([lbl([1 end],:) lbl(:,[1 end])'])));
%     rp = regionprops(holes);
%     min_hole_area = min([rp.Area]);
%     max_hole_area = max([rp.Area]);
% 
%     imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    mask = filteredFaceMask;
% figure
% filteredFaceMask = ~filteredFaceMask;
% im = face;
% s = regionprops(filteredFaceMask, 'Area', 'PixelList');
% [~,ind] = max([s.Area]);
% pix = sub2ind(size(im), s(ind).PixelList(:,2), s(ind).PixelList(:,1));
% out = zeros(size(im));
% out(pix) = im(pix);
% imshow(out);
% % out = imfill(out, 'holes');
% title('funky');
% pause;
%     title('filteredFaceMask');
%     imshow(filteredFaceMask);
%     pause;
%     
%     filteredFaceMask = repmat(filteredFaceMask, [1,1,3]);
%     filteredFaceMask = face.*uint8(filteredFaceMask);
%     
%      title('filteredFaceMask');
%     imshow(filteredFaceMask);
%     pause;
%     
%     filteredFaceMask = imfilter(im2bw(face), fspecial('log'));
% %     filteredFaceMask = imerode(filteredFaceMask, strel('disk',1));
%     filteredFaceMask = imdilate(filteredFaceMask, strel('disk',1));
% %     filteredFaceMask = imerode(filteredFaceMask, strel('disk',2));
%         filteredFaceMask = imfill(filteredFaceMask, 'holes');
%     title('filteredFaceMask');
%     imshow(filteredFaceMask);
%     pause;
    
%     faceHSV = rgb2hsv(face);
%     faceHSV = rgb2ycbcr(face);
    faceHSV = face;
    
%     figure;
%     title('faceHSV');
%     imshow(faceHSV);



% INNAN:
%     dominantColor = [round(mean2(faceHSV(:,:,1))) ...
%                      round(mean2(faceHSV(:,:,2))) ...
%                      round(mean2(faceHSV(:,:,3)))];
%                
%     redChannel = faceHSV(:, :, 1);
%     greenChannel = faceHSV(:, :, 2);
%     blueChannel = faceHSV(:, :, 3);
%     differnece = 1;
%     mask = (abs(redChannel - dominantColor(1,1)) < differnece);
%     mask = or(mask, abs(greenChannel - dominantColor(1,2)) < differnece);
%     mask = or(mask, (abs(blueChannel - dominantColor(1,3)) < differnece));
%     
%     mask = ~mask;



%     figure;
%     title('mask');
%     imshow(mask); 
%     pause;

    
    
    erodeKernel = strel('disk', 5);

    dilateKernel = strel('disk', 15); 
    % maskR = imdilate(maskR, dilateKernel);
   
    
%     mask = imerode(mask, erodeKernel);
%         mask = ExtractNLargestBlobs(mask, 1);
%     mask = imdilate(mask, dilateKernel);
%      figure;
%     title('mask');
%     imshow(mask);
    
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
    dominantHsv = [round(255*mean2(faceHSV(:,:,1))) ...
                   round(255*mean2(faceHSV(:,:,2))) ...
                   round(255*mean2(faceHSV(:,:,3)))];
%                dominantHsv = dominantHsv ./ dominantHsv(1,3);
    dominantHsv = (dominantHsv - min(dominantHsv))./(max(dominantHsv) - min(dominantHsv));
    figure
    imshow(faceHSV);
    title('faceHSV');
%     
%     filtered = imfilter(faceHSV, fspecial('log'));
%     figure
%     imshow(filtered);
%     title('filtered');
%     pause;

% filtered = imfilter(faceHSV, fspecial('sobel'));
h = faceHSV(:,:,1);
s = faceHSV(:,:,2);
v = faceHSV(:,:,3);
% h = h ./ v;
% s = s ./ v;
% 
% minna = [ max(h(:)), min(s(:)), min(v(:)) ];

% max(h(:))
% maskR = h > max(h(:))*0.7;
% maskR = and(maskR, s > max(s(:))*0.7);
% erodeKernel = strel('disk', 1);
% dilateKernel = strel('disk', 5); 
% % maskR = imdilate(maskR, dilateKernel);
% % maskR = imerode(maskR, erodeKernel);
% maskR = imdilate(maskR, dilateKernel);
% maskR = ExtractNLargestBlobs(maskR, 2);
% figure
% imagesc(maskR); 
% colormap gray; 
% % imshow(maskR);
% title('maskR');

% 
% maskR2 = h > max(h(:))*0.99;
% maskR2 = ExtractNLargestBlobs(maskR2, 1);
% figure
% imagesc(maskR2); 
% colormap gray; 
% % imshow(maskR);
% title('maskR2');



% minnas = min(faceHSV(:))
          

%    maskedRgbImage = cat(3, maskedRed, maskedGreen, maskedBlue);
%     rgbImg = and(faceMask, input);
%     rgbImg(:,:,1) = dominantColor;
%     rgbImg(:,:,2) = green;
%     rgbImg(:,:,3) = blue;
    
%     figure
%     imshow(maskedRgbImage);
%     title('maskedRgbImage');
%     pause;
    
%         
%     figure;
%     title('faceMask');
%     imshow(faceMask);
   
     
%     figure;
%     imshow(face);
%     title('face');
    
    [row, col] = find(face(:,:,1) ~= 0);
   
    minCol = min(col);
    maxCol = max(col);
    
    minRow = min(row);
    maxRow = max(row);
    
    cropImg = input(minRow:maxRow, minCol:maxCol, :);
<<<<<<< HEAD
    faceCropped = face(minRow:maxRow, minCol:maxCol, :);

%     % EyeMap & MouthMap
%     faceYCbCr = rgb2ycbcr(face);
%     
% %     figure;
% %     imshow(faceYCbCr);
% %     title('faceYCbCr');
%      
%     faceY = double(faceYCbCr(:,:,1));
%     faceCB = double(faceYCbCr(:,:,2));
%     faceCR = double(faceYCbCr(:,:,3)); 
% 
%     sqFaceCB = faceCB.^2;
%     sqFaceCR = faceCR.^2;
%     cHat = 255 - faceCR;
%     sqCHat = cHat.^2;
%     divFace = faceCB./faceCR;
%     
%     % Now it's normalized to [0,1]
%     % If we need [0,255], just multiply with 255.
%     normFaceCB = (sqFaceCB)/(max(sqFaceCB(:)));
%     normFaceCR = (sqFaceCR)/(max(sqFaceCR(:)));
%     
%     normSqCHat = (sqCHat)/(max(sqCHat(:)));
%     normDivFace = (divFace)/(max(divFace(:)));
%     
%     eyeMapC = 1/3*(normFaceCB + normSqCHat + normDivFace);
%     
% %     figure;
% %     imshow(eyeMapC);
% %     title('eyeMapC');
%     
%     eyeMapCEnhanc = imadjust(eyeMapC);
% 
% %     pause;
% %     erodeKernel = strel('disk', 5);
%     eyeMapCEnhanc = eyeMapCEnhanc - xor(eyeMapCEnhanc, faceMask);
% %     dilateKernel = strel('disk', 10);    
%     
%     eyeMapCEnhanc = eyeMapCEnhanc > 0.9 * max(eyeMapC(:));
    
%     figure;
%     imshow(eyeMapCEnhanc);
%     title('eyeMapCEnhanc');
    

    
    i=face;
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
    imshow(J)
    title('Equalized image');
%     figure
%     J = imadjust(J);
%     J = J > 0.8;
%     imshow(J); title('Equalized image');
%     figure
    
    SE=strel('disk',15,8);
    o=imdilate(igray,SE);
    p=1+imerode(igray,SE);
    eyemaplum=o./p;
    eyemaplum = imadjust(eyemaplum);
    eyemapchr = imadjust(eyemapchr);
    subplot(4,4,9)
    imshow(eyemaplum)
    title('Lum Eyemap');

%     figure
%     imshow(eyemapchr); title('eyemapchr Eyemap'); pause;
%     imshow(eyemapchr); title('eyemaplum Eyemap'); pause;
    
    eyemapchr = eyemapchr > 0.9;
    eyemaplum = eyemaplum > 0.95;
%     figure
%         imshow(eyemapchr); title('eyemapchr Eyemap'); pause;
%     imshow(eyemapchr); title('eyemaplum Eyemap'); pause;
    
    eyemapchr = and(eyemapchr, faceMask);
    eyemaplum = and(eyemaplum, faceMask);
    
    
    cr2 = (cr).^2;
    n = length(cr2(:));
    ovan = (1 / n) * sum(cr2(:));
    cddcb = cr./cb;
    nedan = (1 / n) * sum(cddcb(:));
    
    figure
    nn = 0.95 * (ovan / nedan)
    mouthMap = cr2.*(cr2 - nn*cddcb).^2;
    mouthMap = imadjust(mouthMap);
%     imshow(mouthMap); title('mouthMap'); pause;
    mouthMap = mouthMap > 0.95;
    mouthMap = imfilter(mouthMap, fspecial('sobel'));
    mouthMap = ExtractNLargestBlobs(mouthMap, 1);
     
%     mouthMap = imerode(mouthMap, strel('disk', 1));
    mouthMap = imdilate(mouthMap, strel('disk', 8));
%     imshow(mouthMap); title('mouthMap'); pause;
%     dilateKernel = strel('disk', 4);
%     erodeKernel = strel('disk', 4);
%     mouthMap = imerode(mouthMap, erodeKernel);
%     imshow(mouthMap); title('mouthMap'); pause;
%     mouthMap = imdilate(mouthMap, dilateKernel);
%     imshow(mouthMap); title('mouthMap'); pause;
%     dilateKernel = strel('disk', 5);
%     mouthMap = imdilate(mouthMap, dilateKernel);
%     imshow(mouthMap); title('mouthMap'); pause;
%     mouthMap = imfill(mouthMap, 'holes');
%     imshow(mouthMap); title('mouthMap'); pause;
%     mouthMap = ExtractNLargestBlobs(mouthMap, 1);
%     imshow(mouthMap); title('mouthMap'); pause;
        
    eyeMap = and(eyemapchr, eyemaplum) - mouthMap;
%     figure
    

    mouthMap = repmat(mouthMap, [1,1,3]);
    mouthMap = face.*uint8(mouthMap);
    
    figure;
    imshow(mouthMap);
    title('mouthMap');
    
    dilateKernel = strel('disk', 14);

    figure
%     imshow(eyeMap); title('final Eyemap'); pause;
    
%      for n = 1:6
%         erodeKernel = strel('disk', n);
%         dilateKernel = strel('disk', 1+n);
%         eyeMap = imerode(eyeMap, erodeKernel);
%         eyeMap = imdilate(eyeMap, dilateKernel);
%         
%         title('eyeMap');
%         imshow(eyeMap);
%         pause;
%     end
%     eyeMap = imerode(eyeMap, strel('disk', 2));
%     imshow(eyeMap); title('eyemap'); pause;
%     eyeMap = imdilate(eyeMap, strel('disk', 5));
    
    eyeMap = imfill(eyeMap, 'holes');
%     imshow(eyeMap); title('eyemap'); pause;
%     eyeMap = imerode(eyeMap, strel('disk', 2));
    eyeMap = ExtractNLargestBlobs(eyeMap, 2);
%     imshow(eyeMap); title('eyemap'); pause;
%     eyeMap = imdilate(eyeMap, strel('disk', 2));
%     imshow(eyeMap); title('eyemap'); pause;
%     eyeMap = imdilate(eyeMap, strel('disk', 10));
%     eyeMap = imerode(eyeMap, strel('disk', 10));
    eyeMap = imdilate(eyeMap, strel('disk', 4));
%     imshow(eyeMap); title('final Eyemap'); pause;
    erodeKernel = strel('disk', 11);
%     imshow(eyeMap); title('final Eyemap'); pause;
    eyeMap = imerode(eyeMap, strel('disk', 4));
%     imshow(eyeMap); title('final Eyemap'); pause;
    eyeMap = imdilate(eyeMap, dilateKernel);
%     imshow(eyeMap); title('final Eyemap'); pause;
    eyeMap = imfill(eyeMap, 'holes');
%     imshow(eyeMap); title('final Eyemap'); pause;

    eyeMap = ExtractNLargestBlobs(eyeMap, 2);
%     imshow(eyeMap); title('final Eyemap'); pause;
    finalMask = eyeMap;
    eyeMap = repmat(eyeMap, [1,1,3]);
    eyeMap = face.*uint8(eyeMap);
    imshow(eyeMap); title('final Eyemap'); 

    I=im2bw(eyeMap);  % convert to gray scale
    Rmin=4; Rmax=40;  % circle radius range
    [centersDark, radiiDark, metric] = imfindcircles(I, [Rmin Rmax], ...
                 'ObjectPolarity','dark','sensitivity',0.99)
             
%     centersStrong5 = centersDark(1:2, :)
%     radiiStrong5 = radiiDark(1:2)
%     metricStrong5 = metric(1:2);
    
%     imagesc(I),hold on
%     viscircles(centersStrong5, radiiStrong5,'EdgeColor','b');
    
    eyeMapSize = size(eyeMap)
    imageSizeX = eyeMapSize(1,1,1)
    imageSizeY = eyeMapSize(1,2,1)
    [cols rows ] = meshgrid(1:imageSizeY, 1:imageSizeX);
%     irisMask = final;
    centerX = centersDark(1,1,1)
    centerY = centersDark(1,2,1)
    radius = radiiDark(1,1)
    irisMask = ((rows - centerY).^2 + (cols - centerX).^2) <= radius.^2;
    
    
%     irisMask2 = imdilate(irisMask, strel('r', [5, 2]));
%     irisMapRep = repmat(irisMask2, [1,1,3]);
%     I2 = face;
%     I2(irisMapRep) = 1;
%     imshow(I2); title('I2 Eyemap'); pause;
%     
%     [centersDark, radiiDark, metric] = imfindcircles(I, [Rmin Rmax], ...
%              'ObjectPolarity','dark','sensitivity',0.99)
         
    centerX = centersDark(2,1,1)
    centerY = centersDark(2,2,1)
    radius = radiiDark(2,1)
    irisMask = irisMask | ((rows - centerY).^2 + (cols - centerX).^2) <= radius.^2;
    size(irisMask)
    
    irisMap = repmat(irisMask, [1,1,3]);
    irisMap = face.*uint8(irisMap);
    
    figure
    imshow(irisMap);
    title('irisMap');
    
    eyeColor = [mean2(nonzeros(irisMap(:,:,1))) ...
                mean2(nonzeros(irisMap(:,:,2))) ...
                mean2(nonzeros(irisMap(:,:,3)))]
   
    eyeColor = rgb2hsv(eyeColor/255);


%     sizes = size(centersDark)
   
%     centersDark = find(faceMaskCopy(round(centersDark(1: sizes(1,1), 1), ...
%            faceMaskCopy(round(centersDark(1: sizes(1,1), 2),1)) ~= 0))
%     viscircles(centersDark, radiiDark,'LineStyle','--');hold off



        
    
    
    eyesMouth = or(eyeMap(:,:,1), mouthMap(:,:,1));
    eyesMouth = repmat(eyesMouth, [1,1,3]);
    eyesMouth = face.*uint8(eyesMouth);
    
    figure;
    imshow(eyesMouth);
    title('eyesMouth');

%    eyeMapCEnhanc= houghcircles(eyeMapCEnhanc, 3, 8); 
%    [centers, radii] = imfindcircles(eyeMapCEnhanc, [180, 200], 'Sensitivity', .99);
%     eyeMapCEnhanc= viscircles(centers, radii);

%     eyeMapCEnhanc = imerode(eyeMapCEnhanc, erodeKernel);
%     eyeMapCEnhanc = imdilate(eyeMapCEnhanc, dilateKernel);
%     eyeMapCEnhanc = ExtractNLargestBlobs(eyeMapCEnhanc, 2);
% 
%     
%     diff = abs(normFaceCB - normDivFace);
%     mouthMap = and(diff, normFaceCR);
    
%     figure;
%     imshow(eyeMapCEnhanc);
%     title('eyeMapCEnhanc');

%     output = face;
   
%     I = rgb2gray(imread('face22.jpg'));
%     I = rgb2gray(faceCropped);
    I = rgb2gray(cropImg);
    sizeI = size(I);
    maxSize = max(sizeI);
%     I = imcrop(I,[1 1 maxSize maxSize]);
      
%     I = edge(I, 'sobel');
%      I = imfilter(I, fspecial('disk'));
    I = imfilter(I, fspecial('laplacian'));
%  I = imfilter(I, fspecial('average'));
%  I = imfilter(I, fspecial('sobel'));
% I = imfilter(I, fspecial('log'));
%  size(faceCropped(:,:,1))
%  size(I)
 I = and(faceCropped(:,:,1), I);

%     I = imfilter(I, fspecial('log'));
    
    figure
    imagesc(I); 
    colormap gray; 
    title('Original image'); 

%     imgSize=size(I); %#img is your image matrix
%     finalSize=500;   
%     padImg=zeros(finalSize);

    % padImg(finalSize/2+(1:imgSize(1))-floor(imgSize(1)/2),...
    %        finalSize/2+(1:imgSize(2))-floor(imgSize(2)/2))=I;
    % I = padImge;
    % imagesc(I); colormap gray; title('Original image'); pause;

    [rows, cols] = size(I);

    % Fourier transform
    F = ifftshift(fft2(I))./rows./cols;

    % Show spectrum (log)
%     figure
%     imagesc(log(abs(F))); 
%     title('Fourier transform (abs log)'); 

    % Grid of FFT coordinates
    [rows, cols] = size(F);
    [ux, uy] = meshgrid(([1:cols]-(fix(cols/2)+1))/(cols-mod(cols,2)), ...
                        ([1:rows]-(fix(rows/2)+1))/(rows-mod(rows,2)));

    % Convert to polar coordinates
    th = atan2(uy,ux);
    r = sqrt(ux.^2 + uy.^2);

    % Convert to polar coordinates
    Fr = F .* r;
%     figure
%     imagesc(abs(Fr)); title('Fourier transform x radius'); 
    rcoords = linspace(0, sqrt(ux(1,1)^2 + uy(1,1)^2), rows);
    thcoords = linspace(0, 2*pi, cols);
    [ri, thi] = meshgrid(rcoords, thcoords);
    [x, y] = pol2cart(thi, ri);
    Fp = interp2(ux, uy, abs(Fr), x, y);
%     figure
%     imagesc(Fp); title('Fourier transform in polar coordinates'); 

    % Sum columns to give 1D projection
    F1D = sum(Fp);
    
%     figure
%     plot(rcoords, F1D);
%     title('Projection onto 1D'); 
%     xlim([0 0.5]);
    
%     dydx = diff(F1D(:))./diff(rcoords(:));
    dydx = polyder(F1D);
%     size(polyder(F1D))
%     size(F1D)
%     size(rcoords)

%     figure
%     plot(rcoords(2:end-1), dydx);
%     title('Derivative of Projection onto 1D'); 
%     xlim([0 0.5]);
    
    
%     F1DS = sqrt(ux(1,1)^2 + uy(1,1)^2)*length(F1D)
%     output = [F1D floor(F1DS-5) dominantHsv(1) dominantHsv(2) dominantHsv(3) eyeColor(1) eyeColor(2) eyeColor(3)];

    F1DS = sqrt(ux(1,1)^2 + uy(1,1)^2)*length(dydx);
    output = [dydx floor(F1DS-5) dominantHsv(1) dominantHsv(2) dominantHsv(3) eyeColor(1) eyeColor(2) eyeColor(3)];
    
end
