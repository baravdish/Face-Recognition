clc
close all

addpath src
% 
[access_images, number_of_access_images] = readAllFromDir('access', 'img/access/', '*.jpg');
[no_access_images, number_of_no_access_images] = readAllFromDir('no_access', 'img/no_access/', '*.jpg');
[hard_images, number_of_hard_images] = readAllFromDir('hard', 'img/hard/', '*.jpg');
[tests_images, number_of_tests_images] = readAllFromDir('tests', 'img/tests/17/', '*.jpg');

% result = tnm034(access_images{1})
% result = tnm034(access_images{2})
% result = tnm034(access_images{3})
% result = tnm034(access_images{4})
% return7

warning off
rng(0,'twister');
% im = hard_images{2};

% im = imread('img/tests/17/article-0-1711312F000005DC-491_306x423.jpg');
% balanced_image = colorCorrection(im);
% face_image = detectFace(balanced_image, im);
% 
% im = imread('img/tests/17/maxresdefault.jpg');
% balanced_image = colorCorrection(im);
% face_image = detectFace(balanced_image, im);
% 
im = imread('img/tests/17/pittchanel_01.jpg');
balanced_image = colorCorrection(im);
face_image = detectFace(balanced_image);
% 
% im = imread('img/tests/17/zach-leatherman.jpg');
% balanced_image = colorCorrection(im);
% face_image = detectFace(balanced_image, im);
% 
% im = imread('img/tests/17/oval-face-shape-hairstyles.jpg');
% balanced_image = colorCorrection(im);
% face_image = detectFace(balanced_image, im);

% im = imread('img/tests/17/unnamed.jpg');
% balanced_image = colorCorrection(im);
% face_image = detectFace(balanced_image, im);

% im = no_access_images{4};
% im = no_access_images{2};
% im = hard_images{21};
% im = hard_images{6};
% imshow(im)
% im = tests_images{7};

% im = imread('img/tests/18/url.jpeg');
% balanced_image = colorCorrection(im);
% face_image = detectFace(balanced_image, im);


% result = face_image;
% result = tnm034(im)

disp 'hello'
return;

%%
% for i = 6:6
for i = 1:length(access_images)
    randNumber = randi([0 5]);
    img = access_images{i};
    imgRot = imrotate(img, randNumber, 'crop');
    
    imgScale = imresize(img, 1.1);
    imgHSV = rgb2hsv(img);
    imgHSV(:,:,3) = imgHSV(:,:,3)*0.7;
    imgTone = hsv2rgb(imgHSV);
    imgTone = uint8(255*imgTone);
%     figure
%     imshow(imgTone);

    try
        result = tnm034(imgTone);
    catch
        warning('Program failed.');
        result = 0;
    end
end
%%
for i = 1:length(hard_images)
    
    result = tnm034(hard_images{i})
    
end

%% 
for i = 1:length(no_access_images)
    
    result = tnm034(no_access_images{i})
    
end
