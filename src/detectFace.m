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
    
%     
%     figure; imshow(faceMask); title('faceMask'); pause;
%     faceMask = imerode(faceMask, strel('disk', 10));
%     figure; imshow(faceMask); title('faceMask'); pause;
%     faceMask = ExtractNLargestBlobs(faceMask, 1);
%     figure; imshow(faceMask); title('faceMask'); pause;
%     faceMask = imdilate(faceMask, strel('disk', 10));
%     figure; imshow(faceMask); title('faceMask'); pause;
%     faceMask = imfill(faceMask, 'holes');
%     figure; imshow(faceMask); title('faceMask'); pause;
%     eyesMask = ~imerode(faceMask, strel('disk', 10));
    
    eyesMask =  Cb > Cr + 10;

    imgHSV = rgb2hsv(rgbImage);
    H = imgHSV(:,:,1);
    S = imgHSV(:,:,2);
    V = imgHSV(:,:,3);
    
    hs = V < S - 0.2;
    
%     hs2 = S > V + 0.5;
%     figure; imshow(hs2); title('hs2'); pause;
     
%     figure; imshow(eyesMask); title('eyesMask'); pause;
    
    noFaceMask =  hs | eyesMask;
%     noFaceMask = imfill(noFaceMask, 'holes');
    noFaceMask = imerode(noFaceMask, strel('disk', 5));
    noFaceMask = bwareaopen(noFaceMask, 1000);
    noFaceMask = imdilate(noFaceMask, strel('disk', 5));
%     noFaceMask = imerode(noFaceMask, strel('disk', 1));
    
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
    
    faceMaskRep = repmat(faceMask, [1,1,3]);
    face = rgbImage.*uint8(faceMaskRep);
%     figure; imshow(face); title('face'); pause;
    
%     
%     skinMask2 = faceMask & meanColorize(face, faceMask, 0.01);
%     figure; imshow(skinMask2); title('skinMask2'); pause;

%     faceMaskRep2 = repmat(faceMask, [1,1,3]);
%     face = rgbImage.*uint8(faceMaskRep2);
     
%     faceHSV = rgb2hsv(face);
%     H = faceHSV(:,:,1);
%     S = faceHSV(:,:,2);
%     V = faceHSV(:,:,3);
% 
%     dominantHsv = [round(255*mean2(faceHSV(:,:,1))) ...
%                    round(255*mean2(faceHSV(:,:,2))) ...
%                    round(255*mean2(faceHSV(:,:,3)))];
%     dominantHsv = (dominantHsv - min(dominantHsv))./(max(dominantHsv) - min(dominantHsv));
%     figure
%     imshow(faceHSV);
%     title('faceHSV');pause;
    

%     level = graythresh(face);
%     gray =   im2bw(face, level);
%     filteredFaceMask = edge(rgb2gray(face));

    filteredFaceMask = imfilter(im2bw(imadjust(rgb2gray(face)), 0.47), fspecial('laplacian'));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
%     filteredFaceMask = imfill(filteredFaceMask, 'holes');
    filteredFaceMask = bwareaopen(filteredFaceMask, 50);
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
%     [labeledImage, numberOfBlobs] = bwlabel(filteredFaceMask);
%     blobMeasurements = regionprops(labeledImage, 'area');
% 
%     allAreas = [blobMeasurements.Area];
%     [sortedAreas, sortIndexes] = sort(allAreas, 'descend');
%     
    [labeledImage, numberOfBlobs] = bwlabel(filteredFaceMask);
    convexProperties = regionprops(labeledImage, 'ConvexArea')
    convexAreas = [convexProperties.ConvexArea];
    [sortedConvexAreas, sortedIndices] = sort(convexAreas, 2, 'descend');
    
    filteredFaceMask = ismember(labeledImage, sortedIndices(1));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask'); pause;
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk', 30));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask4'); pause;
    filteredFaceMask2 = filteredFaceMask;
    filteredFaceMask2 = imfill(filteredFaceMask2, 'holes');
    filteredFaceMask = ~xor(~filteredFaceMask, filteredFaceMask2);
    
%     figure; imshow(filteredFaceMask); title('filteredFaceMask5'); pause;

     
%     [labeledMask1 numberOfBlobs1] = bwlabel(filteredFaceMask);
%     measurements1  = regionprops(labeledMask1,'centroid');
%     centroids1 = cat(1, measurements1.Centroid)
%     

%      filteredFaceMask = and(filteredFaceMask, faceMask);
%      
%      [labeledMask2 numberOfBlobs2] = bwlabel(filteredFaceMask);
%      measurements2  = regionprops(labeledMask2,'centroid', 'area');
%      centroids2 = cat(1, measurements2.Centroid)
%      areas2 = [measurements2.Area];
% 
%      indices = []
%     for n=1:length(centroids2)
%         indices(n) = n;
%     end
%     A = areas2;
%     B = indices;
%     A
%     B
%     [sortedA, ind] = sort(A, 2, 'descend');
%     
%     B = ind
%     sortedA
%  
%     
%     newCentroids2 = centroids2
%     for r = 1:size(centroids2, 1)
%         newCentroids2(r, 1) = centroids2(B(r), 1);
%         newCentroids2(r, 2) = centroids2(B(r), 2);
%     end
%     centroids2 = newCentroids2
%     
%      for n2=1:size(centroids2, 1)
%          x2 = centroids2(n2,1);
%          y2 = centroids2(n2,2);
%          for n1=1:size(centroids1, 1)
%             x1 = centroids1(n1,1);
%             y1 = centroids1(n1,2);
%             if (x2 == x1) && (y2 == y1)
%                 filteredFaceMask = ismember(labeledMask2, n2);
%                 n1 = size(centroids1, 1)+1
%                 n2 = size(centroids2, 1)+1;
%                 break
%             end
%          end
%      end

%     figure; imshow(filteredFaceMask); title('filteredFaceMask6'); pause;
    
    
    filteredFaceMask = ExtractNLargestBlobs(filteredFaceMask, 1);
    
    
%     figure; imshow(filteredFaceMask); title('filteredFaceMask7'); pause;
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask9'); pause;
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask10'); pause;
    filteredFaceMask = imdilate(filteredFaceMask, strel('disk',20));
%     figure; imshow(filteredFaceMask); title('filteredFaceMask11'); pause;
    
    faceMask = and(filteredFaceMask, faceMask(:,:,1));
    faceMaskRep = repmat(faceMask, [1,1,3]);
    face = rgbImage.*uint8(faceMaskRep);

    figure; imshow(face); title('face'); pause;
     
%     faceMask = imfill(faceMask, 'holes');
%     faceMask = ExtractNLargestBlobs(faceMask, 1);
%     faceMask = imdilate(faceMask, strel('disk', 20));
%     faceMask = imerode(faceMask, strel('disk', 25));
%     faceMask = ExtractNLargestBlobs(faceMask, 1);
%     
%     skinMask = faceMask & meanColorize(rgbImage, faceMask, 10);
%     figure; imshow(skinMask); title('skinMask'); pause;
%     
%     faceMask = and(faceMask, skinMask);
%     faceMask = imfill(faceMask, 'holes');
%      faceMask = ExtractNLargestBlobs(faceMask, 1);
%     figure; imshow(faceMask); title('faceMask'); pause;

%     faceGradient = imgradient(grayImage) / 255;
%     figure; imshow(faceGradient); title('faceGradient'); pause;
%     faceGradientMasked = faceGradient.*faceMask;
%     figure; imshow(faceGradientMasked); title('faceGradientMasked'); pause;
    
%     newMask = ExtractNLargestBlobs(skinMask, 1);
%     newMask = imdilate(newMask, strel('disk', 4));
%     newMask = imfill(newMask, 'holes');
%     newMask = imerode(newMask, strel('disk', 7));
%     
%     faceMask = and(faceMask, newMask);
%      
%     figure; imshow(faceMask); title('faceMask'); pause;
% 
%     i=rgbImage;
%     figure;
%     subplot(4,4,1)
%     imshow(i)
%     title('original image');
%     iycbcr=rgb2ycbcr(i);
%     iycbcr = im2double(iycbcr);
%     subplot(4,4,2)
%     imshow(iycbcr)
%     title('YCBCR space');
%     y=iycbcr(:,:,1);
%     cb=iycbcr(:,:,2);
%     cr=iycbcr(:,:,3);
%     ccb=cb.^2;
%     subplot(4,4,3)
%     imshow(ccb)
%     title('CB^2');
%     ccr=(1-cr).^2;
%     subplot(4,4,4)
%     imshow(ccr)
%     title('(1-CR)^2');
%     cbcr=ccb./cr;
%     subplot(4,4,5)
%     imshow(cbcr)
%     title('CB/CR');
%     igray=rgb2gray(i);
%     subplot(4,4,6)
%     imshow(igray)
%     title('Gray space');
%     g=1./3;
%     l=g*ccb;
%     m=g*ccr;
%     n=g*cbcr;
% 
%     eyemapchr=l+m+n;
%     subplot(4,4,7)
% 
%     imshow(eyemapchr)
%     title('Chrom Eyemap');
%     J=histeq(eyemapchr);
%     subplot(4,4,8)
%     eyemapchr = J;
% 
%     imshow(J)
%     title('Equalized image');
%     
%     SE=strel('disk',15,8);
%     o=imdilate(igray,SE);
%     p=1+imerode(igray,SE);
%     eyemaplum=o./p;
%     
% 
%     subplot(4,4,9)
%     imshow(eyemaplum)
%     title('Lum Eyemap');
%     pause;
%     
%     figure
%     imshow(eyemapchr); title('eyemapchr Eyemap'); pause;
%     imshow(eyemaplum); title('eyemaplum Eyemap'); pause;



    
    
    output = rgbImage;
end
