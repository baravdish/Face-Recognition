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
    
    figure
    imshow(input);
    title('input');
%     pause;
%     size(input)
    
%     input = create_padded_image(input, 0);
%     input = padarray(input,[16 16],0,'both');
    
%     size(input)
%     title('input');
%     imshow(input);
%     pause;
%     
    imgYCbCr = rgb2ycbcr(input);
    
    Y = imgYCbCr(:,:,1);
    CB = imgYCbCr(:,:,2);
    CR = imgYCbCr(:,:,3);
% 
%     faceMask = and(and(CB > 105, CB < 145), ... 
%                    and(CR > 136, CR < 165));
               
%     faceMask = and(and(CB > 100, CB < 145), ... 
%                    and(CR > 136, CR < 165));

%   ycbcr
% https://web.stanford.edu/class/ee368/Project_03/Project/reports/ee368group15.pdf
%     faceMask = and(and(CB > 100, CB < 145), ... 
%                    and(CR > 132, CR < 165));    
    faceMask = and(and(CB > 100, CB < 145), ... 
                   and(CR > 132, CR < 165));
          
    originalFaceMask = faceMask;
    
    foregroundMask = detractBackground(input);
    
    faceMask = and(faceMask, foregroundMask);
       
    imshow(faceMask); title('faceMask'); pause;
%     faceMask = ExtractNLargestBlobs(faceMask, 1);
%     imshow(faceMask); title('faceMask'); pause;
%     faceMask = meanColorize(input, faceMask, 12);
%     figure; imshow(faceMask); title('meanColorize faceMask'); pause;

%      faceMask = imdilate(faceMask, strel('disk',10));
%       faceMask = imfill(faceMask, 'holes');
%       faceMask = imerode(faceMask, strel('disk',10));
%    figure;  imshow(faceMask); title('faceMask'); pause;
%     faceMask = imdilate(faceMask, strel('disk',2));
%     faceMask = imfill(faceMask, 'holes');
%     faceMask = imerode(faceMask, strel('disk',20));
%     faceMask = imdilate(faceMask, strel('disk',20));
%     imshow(faceMask); title('faceMask'); pause;
               
%     rgb
%     I = input;
%     if(size(I, 3) > 1)
%         final_image = zeros(size(I,1), size(I,2));
%         for i = 1:size(I,1)
%             for j = 1:size(I,2)
%                 R = I(i,j,1);
%                 G = I(i,j,2);
%                 B = I(i,j,3);
%                 if(R > 95 && G > 40 && B > 20)
%                     v = [R,G,B];
%                     if((max(v) - min(v)) > 15)
%                         if(abs(R-G) > 15 && R > G && R > B)
%                             %it is a skin
%                             final_image(i,j) = 1;
%                         end
%                     end
%                 end
%             end
%         end
%     end 
%     figure; imshow(final_image); title('RGB final_image'); pause;
%     faceMask = and(final_image, faceMask);
%     figure; imshow(faceMask); title('and  final_image, faceMask'); pause;
%     
%         faceMask = ExtractNLargestBlobs(faceMask, 1);
%     figure; imshow(faceMask); title('largest blob'); pause;
%     
%     faceMask = imdilate(faceMask, strel('disk',10));
%     figure; imshow(faceMask); title('imdilate blob'); pause;
%     faceMask = imfill(faceMask, 'holes');
%     figure; imshow(faceMask); title('imfill blob'); pause;
%     faceMask = imerode(faceMask, strel('disk',20));
%     figure; imshow(faceMask); title('imerode blob'); pause;
%     faceMask = ExtractNLargestBlobs(faceMask, 1);
%     figure; imshow(faceMask); title('ExtractNLargestBlobs blob'); pause;
%     faceMask = imdilate(faceMask, strel('disk',25));
%     figure; imshow(faceMask); title('imdilate blob'); pause;
%     
%     figure; imshow(faceMask); title('and  final_image, faceMask'); pause;
%     bw = imerode(faceMask, strel('disk', 1));
%     figure; imshow(bw); title('erode and  final_image, faceMask'); pause;
%     bw = bwareaopen(faceMask, 50);
%     figure; imshow(bw); title('open and  final_image, faceMask'); pause;
%     faceMask = meanColorize(input, bw, 10);
%     figure; imshow(faceMask); title('mean again'); pause;
% I=input;
% I=double(I);
% [hue,s,v] = rgb2hsv(I);
% 
% cb =  0.148* I(:,:,1) - 0.291* I(:,:,2) + 0.439 * I(:,:,3) + 128;
% cr =  0.439 * I(:,:,1) - 0.368 * I(:,:,2) -0.071 * I(:,:,3) + 128;
% [w h]=size(I(:,:,1));
% 
% for i=1:w
%     for j=1:h            
%         if  140<=cr(i,j) && cr(i,j)<=165 && 140<=cb(i,j) && cb(i,j)<=195 && 0.01<=hue(i,j) && hue(i,j)<=0.1     
%             segment(i,j)=1;            
%         else       
%             segment(i,j)=0;    
%         end    
%     end
% end
% 
% im(:,:,1)=I(:,:,1).*segment;   
% im(:,:,2)=I(:,:,2).*segment; 
% im(:,:,3)=I(:,:,3).*segment; 
% segment = imdilate(segment, strel('disk', 2));
% segment = imfill(segment, 'holes');
% figure,imshow(segment);
% pause;
% faceMask = segment;




%     img_orig=input;
%    img=img_orig; %copy of original image
%    hsv=rgb2hsv(img);
%    h=hsv(:,:,1);
%    s=hsv(:,:,2);
%   
%    [r c v]=find(h>0.25 | s<=0.15 | s>0.9); %non skin
%    numid=size(r,1);
%   
%    for i=1:numid
%        img(r(i),c(i),:)=0;
%    end
%   
%    figure
%    imshow(img);
%   
%   
%    ycbcr segmentation
%    img_ycbcr=img;  %image from the previous segmentation
%    ycbcr=rgb2ycbcr(img_ycbcr);
%    cb=ycbcr(:,:,2);
%    cr=ycbcr(:,:,3);
%   
%    
%     Detect Skin
%     [r,c,v] = find(cb>=77 & cb<=127 & cr>=133 & cr<=173);
%     [r c v] = find(cb<=77 | cb >=127 | cr<=133 | cr>=173);
%     numid = size(r,1);
%    
%     Mark Skin Pixels
%     for i=1:numid
%         img_ycbcr(r(i),c(i),:) = 0;
%        bin(r(i),c(i)) = 1;
%     end
%    
%     figure
%     imshow(img_ycbcr);
%     title('ycbcr segmentation');
%     pause;
%   rgb segmentation
% 
% img_rgb=img_ycbcr;
% r=img_rgb(:,:,1);
% g=img_rgb(:,:,2);
% b=img_rgb(:,:,3);
% 
% 
% [row col v]= find(b>0.79*g-67 & b<0.78*g+42 & b>0.836*g-14 & b<0.836*g+44 ); %non skin pixels
% numid=size(row,1);
% 
% for i=1:numid
%     img_rgb(row(i),col(i),:)=0;
% end
% 
% figure;
% imshow(img_ycbcr);
% title('img_ycbcr');
% figure;
% imshow(img_rgb);
% title('img_rgb');
% pause;


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
%     max(lapImg(:))
%     min(lapImg(:))
    rega = lapImg.*uint8(bw);
    inva = lapImg.*uint8(~bw);
%     max(rega(:))
%     min(inva(:))
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
    
    figure
    imshow(foreground);
    title('foreground');
    pause;
    
  
    
    
               
    figure
    title('faceMask');
    imshow(faceMask);
%     pause;
    faceMask = imfill(faceMask, 'holes');
    
%     title('faceMask');
%     imshow(facsk);
%     pause;
    faceMask = ExtractNLargestBlobs(faceMask, 1);
%     
%     faceMaskRep = reeMapmat(faceMask, [1,1,3]);
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
    
    if homogeneousBackground == 1
          faceMask = and(foreground, faceMask);
    end
    
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
    
    figure
    imshow(faceMask);
    title('cropped faceMask');
    pause;
%     
    figure; imshow(faceMask); title('faceMask');  pause;
    
    newMask = meanColorize(input, faceMask);
    originalNewMask = newMask;
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
    
    figure;
    imshow(face);
    title('face');
    pause;

    
    i=input;
    iMask = face > 0;
    iMask = or(or( iMask(:,:,1), iMask(:,:,2)), iMask(:,:,3));
    iMaskRep = repmat(iMask, [1,1,3]);
%     face = input.*uint8(iMask);
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
%     figure;
    imshow(eyemapchr)
    title('Chrom Eyemap');
    J=histeq(eyemapchr);
    subplot(4,4,8)
    eyemapchr = J;
%     figure
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
    
%     size(eyemaplum)
%     size(iMask)
%     size(iMaskRep)
    
%     figure
%     eyeMap = and(eyemapchr, eyemaplum);
%     imshow(eyeMap); title('eyeMap'); pause;
    
%     eyemaplum = eyemaplum.*uint8(iMask);
%     eyemapchr = eyemapchr.*iMask;
%     eyemaplum = imadjust(eyemaplum);
%     eyemapchr = imadjust(eyemapchr);
    subplot(4,4,9)
%     figure
    imshow(eyemaplum)
    title('Lum Eyemap');
    pause;
    
    
    figure
    imshow(eyemapchr); title('eyemapchr Eyemap'); pause;
    imshow(eyemaplum); title('eyemaplum Eyemap'); pause;
    
    
    
    
    
%     figure;
%     img = eyemaplum;
%     overSaturated = rgb2gray(img);
%     %     imshow(overSaturated); title('Over saturation'); pause;
%     %     imshow(imadjust(overSaturated)); title('imadjust(overSaturated)'); pause;
% %     imshow(imfilter(imadjust(overSaturated), fspecial('laplacian'))); title('filter imadjust(overSaturated)'); pause;
%     overSaturated = imfilter(overSaturated, fspecial('laplacian'));
%     %     imshow(overSaturated); title('filter Over saturation'); pause;
%     overSaturated = imdilate(overSaturated, strel('disk', 2));
% %     imshow(overSaturated); title('dilate Over saturation'); pause;
%     overSaturated = ~overSaturated;
%     imshow(overSaturated); title('filter Over saturation'); pause;
    
%     imshow(img); title('before Over saturation'); pause;
%     overSaturated = ~overSaturated;
%     overSaturatedRep = repmat(overSaturated, [1,1,3]);
%     img = img.*uint8(overSaturatedRep);
%     imshow(img); title('after Over saturation'); pause;
    
%     originalImg = img;

%     subplot(1,4, 1);
%     imshow(originalImg);
%     subplot(1,4, 2);
%     bar(pixelCounts);
%     subplot(1,4, 3);
%     plot(cdf);
%     subplot(1,4, 4);
%     thresholdIndex = find(cdf < 0.96, 1, 'last');
%     thresholdValue = grayLevels(thresholdIndex);
%     thresholdedImage = originalImg;
%     thresholdedImage(originalImg(:) < thresholdValue) = 0;
%     imshow(thresholdedImage);
%     set(gcf, 'Position', get(0, 'ScreenSize')); % Maximize figure.
%     pause;
    
%     figure; imshow(eyemapchr); title('eyemapchr'); pause;
%     [pixelCounts grayLevels] = imhist(eyemapchr);
%     cdf = cumsum(pixelCounts) / sum(pixelCounts);
%     
%     figure; imshow(eyemaplum); title('eyemaplum'); pause;
%     [pixelCounts grayLevels] = imhist(eyemaplum);
%     cdf = cumsum(pixelCounts) / sum(pixelCounts);
%     eyemaplumIndex = find(cdf < 0.99, 1, 'last');
%     eyemaplumTh = grayLevels(eyemaplumIndex) / 255;
    
    vea = eyemapchr.*finalFaceMask;
    imshow(vea); title('vea'); pause;
    vec = nonzeros(vea(:));
    soretedVec = sort(vec, 'descend');
    eyemapchrTh = 0.95 * soretedVec(1,1)
%     soretedVec = sort(vec, 'ascend');
%     len = length(soretedVec);
%     index = 0.95 * len;
%     eyemapchrTh = soretedVec(floor(index));
    
    eyemapchr22 = imgradient(rgb2gray(input)) / 255;
    imshow(eyemapchr22);
    title('imgradient');
    pause;


    figure; imshow(eyemapchr); title('eyemapchr before'); pause;
%     eyemapchr2 = medfilt2(eyemapchr);
%     eyemapchr2 = medfilt2(eyemapchr2);
%     eyemapchr2 = medfilt2(eyemapchr2);
%     eyemapchr = eyemapchr + eyemapchr2;
    %     eyemapchr = eyemapchr.*eyemapchr;
%     eyemapchr = imadjust(eyemapchr);
%     eyemapchr = imdilate(eyemapchr, strel('disk',5));
%     eyemapchr = imfilter(eyemapchr, fspecial('gaussian', 25));
%     level = graythresh(eyemapchr);
%     eyemapchr = im2bw(eyemapchr, 2*level);
%     eyemapchr = bradley(eyemapchr, [125 125], 10);
%     eyemapchr = imfilter(eyemapchr, fspecial('gaussian'));
%     eyemapchr = imfilter(eyemapchr, fspecial('gaussian'));
    figure; imshow(eyemapchr); title('eyemapchr after'); pause;
    eyemapchr = eyemapchr > eyemapchrTh;
%     eyemapchr = and(eyemapchr, ~originalNewMask);
    figure; imshow(eyemapchr); title('eyemapchr'); pause;
    
    
%     size(eyemaplum)
%     size(eyemaplum)
    vea = eyemaplum.*uint8(finalFaceMask);
    imshow(vea); title('vea'); pause;
    vec = nonzeros(vea(:));
    soretedVec = sort(vec, 'descend');
    eyemaplumTh = 0.9 * (soretedVec(1,1) / 255)
%     soretedVec = sort(vec, 'ascend');
%     len = length(soretedVec);
%     index = 0.9 * len;
%     eyemaplumTh = soretedVec(floor(index));
    
    eyemaplum = eyemaplum > eyemaplumTh;
    figure; imshow(eyemaplum); title('eyemaplum'); pause;
% nPixels = length(rVec);
% percentage = 0.001;
% nFind = round(nPixels*percentage);
% meanValues = [mean(sortedImg((1:nFind),1)), ...
%               mean(sortedImg((1:nFind),2)), ...
%               mean(sortedImg((1:nFind),3))];
          
%     figure
%         imshow(eyemapchr); title('eyemapchr Eyemap'); pause;
%     imshow(eyemapchr); title('eyemaplum Eyemap'); pause;
    
    eyemapchr = and(eyemapchr, faceMask);
  
    
    eyemaplum = and(eyemaplum, faceMask);
    
    

    imshow(eyemapchr); title('eyemapchr true'); pause;
    imshow(eyemaplum); title('eyemaplum true'); pause;
    eyeMap = and(eyemapchr, eyemaplum);
    imshow(eyeMap); title('eyeMap'); pause;
%     figure
%     figure; imshow(eyeMap); title('eyeMap1'); pause;
%     figure; imshow(finalFaceMask); title('finalFaceMask'); pause;
    figure;
    overSaturated = rgb2gray(face);
%     imshow(overSaturated); title('Over saturation'); pause;
%     imshow(imadjust(overSaturated)); title('imadjust(overSaturated)'); pause;
    imshow(imfilter(imadjust(overSaturated), fspecial('laplacian'))); title('filter imadjust(overSaturated)'); pause;
    overSaturated = imfilter(overSaturated, fspecial('laplacian'));
%     imshow(overSaturated); title('filter Over saturation'); pause;
    overSaturated = imdilate(overSaturated, strel('disk', 2));
    imshow(overSaturated); title('dilate Over saturation'); pause;
    overSaturated = ~overSaturated;
    imshow(overSaturated); title('filter Over saturation'); pause;
    
    overSaturated = and(finalFaceMask, overSaturated);
    imshow(overSaturated); title('Over saturation'); pause;
    
    imshow(eyeMap); title('eyeMap11'); pause;
    eyeMap = eyeMap - and(overSaturated,eyeMap);
    
    imshow(eyeMap); title('eyeMap2'); pause;
    
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
    eyeMap = imdilate(eyeMap, strel('disk', 4));
%     eyeMap2 = ExtractNLargestBlobs(eyeMap, 2);
%     eyeMap2 = imdilate(eyeMap2, strel('disk', 2));
%     eyeMap = or(eyeMap, eyeMap);
    eyeMap = ExtractNLargestBlobs(eyeMap, 2);

%     imshow(eyeMap); title('eyemap'); pause;
%     eyeMap = imdilate(eyeMap, strel('disk', 2));
%     imshow(eyeMap); title('eyemap'); pause;
%     eyeMap = imdilate(eyeMap, strel('disk', 10));
%     eyeMap = imerode(eyeMap, strel('disk', 10));
    
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
%     imshow(eyeMap); title('final Eyemap'); pause;
    finalEyeMask = eyeMap;
    eyeMap = repmat(eyeMap, [1,1,3]);
    eyeMap = face.*uint8(eyeMap);
    imshow(eyeMap); title('final Eyemap'); 

%     I=im2bw(eyeMap);  % convert to gray scale
    level = graythresh(eyeMap);
    I = im2bw(eyeMap, level);
    imshow(I); title('I Eyemap'); pause;
%     I=rgb2gray(eyeMap);  % convert to gray scale
%     imshow(I); title('I Eyemap'); pause;
%     I=imadjust(rgb2gray(eyeMap));  % convert to gray scale
%     imshow(I); title('I Eyemap'); pause;
    Rmin=4; Rmax=40;  % circle radius range
    [centersDark, radiiDark, metric] = imfindcircles(I, [Rmin Rmax], ...
                 'ObjectPolarity','dark','sensitivity',0.99);
     
    figure
    imshow(eyeMap); title('final Eyemap'); pause;
%     viscircles(centersBright, radiiBright,'EdgeColor','b');
    viscircles(centersDark, radiiDark,'LineStyle','--');
    pause
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
    
    figure
    imshow(irisMap);
    title('irisMap');
    
    eyeColor = [mean2(nonzeros(irisMap(:,:,1))) ...
                mean2(nonzeros(irisMap(:,:,2))) ...
                mean2(nonzeros(irisMap(:,:,3)))]
   
    eyeColor = rgb2hsv(eyeColor/255);

%     figure
%     imshow(face);
%     title('original');
%     pause;
    referenceDirection = [1 0];
    eyeDirection = [centerX1-centerX2, centerY1-centerY2];
    eyeDirection = eyeDirection / norm(eyeDirection);
    angle = acos(dot(referenceDirection, eyeDirection)) * 180 / pi
    angle = min(angle, 180-angle);
    

    figure
    imshow(face);
    title('face before rotated');
    pause;
    face = imrotate(face, -angle);
    eyeMap = imrotate(eyeMap, -angle);
    irisMap = imrotate(irisMap, -angle);
    cropImg = imrotate(cropImg, -angle);
    faceCropped = imrotate(faceCropped, -angle);
    finalFaceMask = imrotate(finalFaceMask, -angle);
    finalEyeMask = imrotate(finalEyeMask, -angle);
    overSaturated = imrotate(overSaturated, -angle);
    figure
    imshow(face);
    title('face after rotated');
    pause;
    
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
    
%     mouthMap = mouthMap.*
    mouthMap = imadjust(mouthMap);
%       mouthMap =  histeq(mouthMap);
%     size(mouthMap)
%     size(finalFaceMask)
%     finalFaceMaskRep = repmat(finalFaceMask, [1,1,3]);
    vea = mouthMap;
    imshow(vea); title('vea'); pause;
    vec = nonzeros(vea(:));
    soretedVec = sort(vec, 'descend');
%     maxa = max(mouthMap(nonzeros(finalFaceMask)))
    mouthTh = 0.99 * soretedVec(1,1)
    pause;
    imshow(mouthMap); title('mouthMap'); pause;
    mouthMap = mouthMap > mouthTh;
%     mouthMap = mouthMap - and(mouthMap, finalEyeMask);
%     mouthMap = mouthMap - and(mouthMap, overSaturated);
%     mouthMap = mouthMap - and(mouthMap, finalFaceMask);
    imshow(mouthMap); title('mouthMap'); pause;
%    mouthMap = edge(mouthMap,'sobel','horizontal');
%     imshow(mouthMap); title('mouthMap'); pause;
    mouthMap = or(imfilter(mouthMap, [-1 0 1]'), imfilter(mouthMap, [1 0 -1]'));
%     mouthMap = or(mouthMap, imfilter(mouthMap, fspecial('sobel')));
%       mouthMap = edge(mouthMap, 'Roberts');
%       mouthMap = edge(mouthMap, 'Canny');
%       'Canny'
        
    mouthMap = uint8(mouthMap).*uint8(finalFaceMask);
    mouthMap = uint8(mouthMap) - uint8(mouthMap).*uint8(finalEyeMask);
    mouthMap = uint8(mouthMap) - uint8(mouthMap).*uint8(overSaturated);
%     mouthMap = mouthMap - and(mouthMap, finalEyeMask);
%     mouthMap = mouthMap - and(mouthMap, overSaturated);
%     imshow(mouthMap); title('mouthMap'); pause;
%     mouthMap = or(imfilter(mouthMap, [-1 0 1]'), ...
%                   imfilter(mouthMap, fspecial('sobel')));
    imshow(mouthMap); title('mouthMap'); pause;
    
%     imshow(mouthMap); title('mouthMap'); pause;
    
%     mouthMap = imdilate(mouthMap, strel('disk', 1));
%     mouthMap2 = ExtractNLargestBlobs(mouthMap, 1);
%     mouthMap2 = imdilate(mouthMap2, strel('disk', 4));
%     mouthMap = or(mouthMap, mouthMap2);
    mouthMap = ExtractNLargestBlobs(mouthMap, 1);
    imshow(mouthMap); title('mouthMap'); pause;
%     mouthMap = imerode(mouthMap, strel('disk', 1));
    mouthMap = imdilate(mouthMap, strel('disk', 8));
    mouthMap = imfill(mouthMap, 'holes');
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

    mouthMap = repmat(mouthMap, [1,1,3]);
    mouthMap = face.*uint8(mouthMap);
    
    figure;
    imshow(mouthMap);
    title('mouthMap');
    
    
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
%     I = imadjust(I);
    sizeI = size(I);
    maxSize = max(sizeI);
%     I = imcrop(I,[1 1 maxSize maxSize]);
I = imfilter(I, fspecial('laplacian'));

%      I = imfilter(I, fspecial('disk', 10));     
%     I = imfilter(I, fspecial('laplacian'));
%     I = imerode(I, strel('disk', 1));
%     I = ExtractNLargestBlobs(I,2);
%     I = imdilate(I, strel('disk', 5));
%     
%     IREP = repmat(I>0, [1,1,3]);
%     Ie = cropImg.*uint8(IREP);
%     
%      figure
%     imshow(Ie); 
%     colormap gray; 
%     title('COOL EYE image'); 
    
    
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
%     COOOOOOOOOOOOOOOOL = size(dydx)
% langd = length(dydx);  
langd = floor(F1DS-5);
val = dydx;
% val = F1D;
output = [val langd dominantHsv(1) dominantHsv(2) dominantHsv(3) eyeColor(1) eyeColor(2) eyeColor(3)];

end
