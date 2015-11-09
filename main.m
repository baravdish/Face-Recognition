%% Normalization part of main.m
warning('off','all')
warning('query','all')

clear
clc
close all

tic
% load DB

nImages = 6;
path = 'DB1\db1_0';
fileformat = '.jpg';

%% Main
% for i = 1:1
for i = 1:nImages
    
    img = imread(strcat(path, num2str(i), fileformat));
    
    imgYCbCr = rgb2ycbcr(img);
    imgHSV = rgb2hsv(img);
    
    Y = imgYCbCr(:,:,1);a
    CB = imgYCbCr(:,:,2);
    CR = imgYCbCr(:,:,3);
    
    H = imgHSV(:,:,1);
    S = imgHSV(:,:,2);
    V = imgHSV(:,:,3);
    
    figure
    fontSize = 12;
    
    subplot(4, 3, 1);
	imshow(Y);
	title('Y', 'FontSize', fontSize);
	subplot(4, 3, 2);
	imshow(CB);
	title('CB', 'FontSize', fontSize);
	subplot(4, 3, 3);
	imshow(CR);
	title('CR', 'FontSize', fontSize);
    
    subplot(4, 3, 4);
	plot(imhist(Y));
	title('hist Y', 'FontSize', fontSize);
	subplot(4, 3, 5);
	plot(imhist(CB));
	title('hist CB', 'FontSize', fontSize);
	subplot(4, 3, 6);
	plot(imhist(CR));
	title('hist CR', 'FontSize', fontSize);
    
    subplot(4, 3, 7);
	imshow(H);
	title('H', 'FontSize', fontSize);
	subplot(4, 3, 8);
	imshow(S);
	title('S', 'FontSize', fontSize);
	subplot(4, 3, 9);
	imshow(V);
	title('V', 'FontSize', fontSize);
    
    subplot(4, 3, 10);
	plot(imhist(H));
	title('hist H', 'FontSize', fontSize);
	subplot(4, 3, 11);
	plot(imhist(S));
	title('hist S', 'FontSize', fontSize);
	subplot(4, 3, 12);
	plot(imhist(V));
	title('hist V', 'FontSize', fontSize);
    
    faceMask = and(and(CB > 105, CB < 135), ...
                   and(CR > 140, CR < 165));
    
%     faceMask = and(H > 6/360, H < 38/360);
    
%     faceMask = and(and(and(CB > 105, CB < 135), ...
%                    and(CR > 140, CR < 165)), ...
%                    and(H > 6/360, H < 39/360));    
                
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

%     figure
%     title('faceMask' + int2str(i))
%     imshow(faceMask);
    
    faceMask = repmat(faceMask, [1,1,3]);
    face = img.*uint8(faceMask);
    
    figure
    title('face' + int2str(i))
    imshow(face);
   
    
end

%%
%     threshCB = and(CB > 105, CB < 135);
%     threshCR = and(CR > 140, CR < 165);
%     
% %   10 erode och 25 dilate verkar vara lämpligt
%     SEerode = strel('square', 10);
%     SEdilate = strel('square', 25);
%     
%     erosionCR = imerode(threshCR, SEerode);
%     closingCR= imdilate(erosionCR, SEdilate);
%     
%     figure(i)
%     imshow(erosionCR);
%     
% %     Matlabs egna closing
% %     closingCR = imclose(closingCR, SEdilate);
%     
%     closeCR3 = repmat(closingCR, [1,1,3]);
%     
%     imgBinar = img;
%     
%     imgBinar(~closeCR3) = 0;
%     
%     binar = imgBinar ~= 0;
%     
%     [row, col] = find(binar(:,:,2) ~= 0);
%     
%     minCol = min(col);
%     maxCol = max(col);
% 
%     minRow = min(row);
%     maxRow = max(row);
%     realImgCrop = img(minRow:maxRow, minCol:maxCol);
% 
%     width = maxRow - minRow;
%     height = maxCol - minCol;
%     croppedImg = imcrop(imgBinar ,[minCol minRow width height]);
    
%     TODO Eyemap:
    
    

%     figure(i)
%     imshow(realImgCrop);


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
