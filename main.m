clc
close all

addpath src

[access_images, number_of_access_images] = readAllFromDir('access', 'img/access/', '*.jpg');
[no_access_images, er_of_no_access_images] = readAllFromDir('no_access', 'img/no_access/', '*.jpg');
[hard_images, number_of_hard_images] = readAllFromDir('hard', 'img/hard/', '*.jpg');

% result = tnm034(access_images{1})
% result = tnm034(access_images{2})
% result = tnm034(access_images{3})
% result = tnm034(access_images{4})
% return

% im = hard_images{15};
% imshow(im)
% result = tnm034(im)
% return;
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
