%% Normalization part of main.m
warning('off','all')
warning('query','all')

clear
clc
close all

tic
% load DB

nImages = 9;
path = 'DB2\bl_0';
fileformat = '.jpg';

%% Main
% for i = 3:3
close all
for i = 1:nImages
    
    img = imread(strcat(path, num2str(i), fileformat));
%     imgGray = rgb2gray(img);
%     imshow(imgGray);
%     hold on
%     h = fspecial('sobel');
%     imgF = imfilter(img, h);
%     figure
%     imshow(I);
%     figure
%     imshow(img)
%     title('img');
    
    imgNorm = lightNorm(img);
%     imgNorm = whiteNorm(img);
    
%     figure
%     imshow(imgNorm);
%     title('imgNorm')
    
    imgYCbCr = rgb2ycbcr(imgNorm);
%     imgHSV = rgb2hsv(img);
    
    Y = imgYCbCr(:,:,1);
    CB = imgYCbCr(:,:,2);
    CR = imgYCbCr(:,:,3);
    
%     H = imgHSV(:,:,1);
%     S = imgHSV(:,:,2);
%     V = imgHSV(:,:,3);
    

%     faceMask = and(and(CB > 105, CB < 135), ... 
%     and(CR > 140, CR < 165)); 

    faceMask = and(and(CB > 100, CB < 140), ... 
                   and(CR > 135, CR < 170)); 


    figure
    imshow(faceMask);
    title('faceMask');

    for n = 1:6
        erodeKernel = strel('disk', n);
        dilateKernel = strel('disk', 1+n);
        faceMask = imerode(faceMask, erodeKernel);
        faceMask = imdilate(faceMask, dilateKernel);
    end
   
    erodeKernel = strel('disk', 25);
    dilateKernel = strel('disk', 20);    
    faceMask = imdilate(faceMask, dilateKernel);
    faceMask = imerode(faceMask, erodeKernel);

    erodeKernel = strel('disk', 10);
    dilateKernel = strel('disk', 10); 
    faceMask = imerode(faceMask, erodeKernel);
    faceMask = imdilate(faceMask, dilateKernel);
    
    faceMask = repmat(faceMask, [1,1,3]);
    face = img.*uint8(faceMask);
    
    figure
    imshow(face);
    title('face')
    
    [row, col] = find(face(:,:,1) ~= 0);
    
    minCol = min(col);
    maxCol = max(col);

    minRow = min(row);
    maxRow = max(row);
    cropImg = img(minRow:maxRow, minCol:maxCol, :);
%     figure
%     imshow(cropImg);
     
    faceYCbCr = rgb2ycbcr(cropImg);
%     
    faceY = double(faceYCbCr(:,:,1));
    faceCB = double(faceYCbCr(:,:,2));
    faceCR = double(faceYCbCr(:,:,3)); 

    sqFaceCB = faceCB.^2;
    sqFaceCR = faceCR.^2;
    cHat = 255 - faceCR;
    sqCHat = cHat.^2;
    divFace = faceCB./faceCR;
    
    % Now it's normalized to [0,1]
    % If we need [0,255], just multiply with 255.
    normFaceCB = (sqFaceCB)/(max(sqFaceCB(:)));
    normFaceCR = (sqFaceCR)/(max(sqFaceCR(:)));
    
    normSqCHat = (sqCHat)/(max(sqCHat(:)));
    normDivFace = (divFace)/(max(divFace(:)));
    
    eyeMapC = 1/3*(normFaceCB + normSqCHat + normDivFace);
    imshow(eyeMapCEnhanc);
    hold on
    eyeMapCEnhanc = imadjust(eyeMapC);
    [centers, radii, metric] = imfindcircles(eyeMapCEnhanc,[1 20]);
    centersStrong5 = centers(1:5,:);
    radiiStrong5 = radii(1:5);
    metricStrong5 = metric(1:5);
    viscircles(centersStrong5, radiiStrong5,'EdgeColor','b');

%     figure
    diff = abs(normFaceCB - normDivFace);
    mouthMap = and(diff, normFaceCR);
%     figure
%     imshow(mouthMap, []);
%        
    
    
%     figure
%     title('face' + int2str(i))
%     imshow(face);
   
    
end

%% Otsu's method
for i = 1:nImages
    
    img = imread(strcat(path, num2str(i), fileformat));
    imgGray = rgb2gray(img);
%     img = album{i};
    normImg = otsuNormalize(imgGray);
        
    normImg3 = repmat(normImg, [1,1,3]);
    
    img(normImg3) = 0;
    
    figure
    imshow(img)

end
% save DB.mat
toc

%% TEST part
% - A test to check the difference between this implementation and the built-in function graythresh()

for n = 1:4
    test = threshImages{n};
    sumMy = sum(test(:));
    a = album{n};
    level = graythresh(a);
    sumReal = a > 255*level;
    sumReal = abs(sum(sumReal(:)));
    difference = abs(sumMy - sumReal);
    difference/(256*256)
end
