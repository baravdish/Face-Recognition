%% Normalization part of main.m
clear
clc

tic
load DB

nImages = 4;
path = 'DB0\db0_';
fileformat = '.jpg';

for i = 1:nImages
    
%     img = imread(strcat(path, num2str(i), fileformat));
%     img = rgb2gray(img);
    img = album{i};
    normImg = otsuNormalize(img);
    figure
    imshow(normImg)
    
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
