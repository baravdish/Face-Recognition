clc
close all

addpath src

[access_images, number_of_access_images] = readAllFromDir('access', 'img/access/', '*.jpg');
[no_access_images, number_of_no_access_images] = readAllFromDir('no_access', 'img/no_access/', '*.jpg');
[hard_images, number_of_hard_images] = readAllFromDir('hard', 'img/hard/', '*.jpg');
[tests_images, number_of_tests_images] = readAllFromDir('tests', 'img/tests/17/', '*.jpg');
warning off

nAccImages = length(access_images);
nNoAccImages = length(no_access_images);
nHardImages = length(hard_images);
nTestImages = length(tests_images);

%%
tic
% En struct som är Nx1 stor med resp. värden
result = repmat(struct('blurFace', 1 , 'blurResult', 1,...
                        'toneFace', 1, 'toneResult', 1, 'scaleFace', 1, ...
                        'scaleResult', 1, 'rotationFace', 1, 'rotationResult', ...
                        1), nHardImages, 1 );

for i = 1:nHardImages
    
    % TODO: Borde ta in alla bilder inte bara access, bunta ihop alla eller
    % en till loop?
%     img = access_images{i};
%     img = no_access_images{i};
    img = hard_images{i};
    
    % PARAMETERS
    sigma = 1.0;
    light = 0.85;
    scale = 1.1;
    rotation = 5;
    height = 200; % Requires that the database is rebuilt!
    width = 200; % Requires that the database is rebuilt!
    threshold = 0.001140; % Requires that the database is rebuilt!
    kernel_size = 21; % Requires that the database is rebuilt!

    imgBlur = applyLP(img, sigma);
    imgBlur = colorCorrection(imgBlur);
    result(i).blurFace = detectFace(imgBlur);
    result(i).blurResult = verify(result(i).blurFace, height, width, threshold, kernel_size, 0, 0);

    imgTone = applyTone(img, light);
    imgTone = colorCorrection(imgTone);
    result(i).toneFace = detectFace(imgTone);
    result(i).toneResult = verify(result(i).toneFace, height, width, threshold, kernel_size, 0, 0);

    imgScale = imresize(img, scale);
    imgScale = colorCorrection(imgScale);
    result(i).scaleFace = detectFace(imgScale);
    result(i).scaleResult = verify(result(i).scaleFace, height, width, threshold, kernel_size, 0, 0);

    imgRot = imrotate(img, rotation, 'crop');
    imgRot = colorCorrection(imgRot);
    result(i).rotationFace = detectFace(imgRot);
    result(i).rotationResult = verify(result(i).rotationFace, height, width, threshold, kernel_size, 0, 0);
    
%     figure; imshow(result(i).blurFace); title(strcat('blurFace, with sigma = ', num2str(sigma)));
%     figure; imshow(result(i).toneFace); title(strcat('toneFace, with tonevalue = ', num2str(light)));
%     figure; imshow(result(i).scaleFace); title(strcat('toneFace, with scale = ', num2str(scale)));
%     figure; imshow(result(i).rotationFace); title(strcat('toneFace, with rotation = ', num2str(rotation), ' degrees'));

end
toc
%% Evaluate matching

blurRes = 0;
toneRes = 0;
scaleRes = 0;
rotationRes = 0;
nImages = length(result);
for i = 1:nImages
    
    if(result(i).blurResult == i)
        blurRes = blurRes + 1;
    elseif(result(i).toneResult == i)
        toneRes = toneRes + 1;
    elseif(result(i).scaleResult == i)
        scaleRes = scaleRes + 1;
    elseif(result(i).scaleResult == i)
        scaleRes = scaleRes + 1;     
    end
    
end


rotationHitRate = rotationRes/nImages;
scaleHitRate = scaleRes/nImages;
blurHitRate = blurRes/nImages;
toneHitRate = toneRes/nImages;






