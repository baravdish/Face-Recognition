%% Normalization part of main.m
clear
clc

tic
% load DB

nImages = 4;
path = 'DB0\db0_';
fileformat = '.jpg';

%% Main


for i = 1:nImages
    
    img = imread(strcat(path, num2str(i), fileformat));
    
    imgCB = rgb2ycbcr(img);
    
    Y = imgCB(:,:,1);
    CB = imgCB(:,:,2);
    CR = imgCB(:,:,3);
    
    threshCB = and(CB > 105, CB < 135);
    threshCR = and(CR > 140, CR < 165);
    
%   10 erode och 25 dilate verkar vara lämpligt
    SEerode = strel('square', 10);
    SEdilate = strel('square', 25);
    
    erosionCR = imerode(threshCR, SEerode);
    closingCR= imdilate(erosionCR, SEdilate);

%     Matlabs egna closing
%     closingCR = imclose(closingCR, SEdilate);
    
    closeCR3 = repmat(closingCR, [1,1,3]);
    
    imgBinar = img;
    
    imgBinar(~closeCR3) = 0;
    
    binar = imgBinar ~= 0;
    
    [row, col] = find(binar(:,:,2) ~= 0);
    
    minCol = min(col);
    maxCol = max(col);

    minRow = min(row);
    maxRow = max(row);
    realImgCrop = img(minRow:maxRow, minCol:maxCol);

    width = maxRow - minRow;
    height = maxCol - minCol;
    croppedImg = imcrop(imgBinar ,[minCol minRow width height]);
    
%     TODO Eyemap:


    figure(i)
    imshow(realImgCrop);
    
    
    
    
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
