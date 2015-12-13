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
                        1, 'sigma', 1, 'light', 1, 'scale', 1, 'rotation', 1), nTestImages, 1 );

height = 200; % Requires that the database is rebuilt!
width = 200; % Requires that the database is rebuilt!
threshold = 0.001140; % Requires that the database is rebuilt!
kernel_size = 21; % Requires that the database is rebuilt!
h = waitbar(0, 'Performing test cases...');
for i = 1:length(result)
    
    % TODO: Borde ta in alla bilder inte bara access, bunta ihop alla eller
    % en till loop?
%     img = access_images{i};
%     img = no_access_images{i};
%     img = hard_images{i};
    img = tests_images{i};
    
    result(i).sigma = 1.0;
    result(i).light = 1.3;
    result(i).scale = 1.1;
    result(i).rotation = 20;

%     imgBlur = applyLP(img, result(i).sigma);
%     imgBlur = colorCorrection(imgBlur);
%     result(i).blurFace = detectFace(imgBlur);
%     result(i).blurResult = verify(result(i).blurFace, height, width, threshold, kernel_size, 0, 0);
% 
%     imgTone = applyTone(img, result(i).light);
%     imgTone = colorCorrection(imgTone);
%     result(i).toneFace = detectFace(imgTone);
%     result(i).toneResult = verify(result(i).toneFace, height, width, threshold, kernel_size, 0, 0);

%     imgScale = imresize(img, result(i).scale);
%     imgScale = colorCorrection(imgScale);
%     result(i).scaleFace = detectFace(imgScale);
%     result(i).scaleResult = verify(result(i).scaleFace, height, width, threshold, kernel_size, 0, 0);

    imgRot = imrotate(img, result(i).rotation, 'crop');
    imgRot = colorCorrection(imgRot);
    result(i).rotationFace = detectFace(imgRot);
%     result(i).rotationResult = verify(result(i).rotationFace, height, width, threshold, kernel_size, 0, 0);
    
%     figure; imshow(result(i).blurFace); title(strcat('blurFace, with sigma = ', num2str(sigma)));
%     figure; imshow(result(i).toneFace); title(strcat('toneFace, with tonevalue = ', num2str(light)));
%     figure; imshow(result(i).scaleFace); title(strcat('toneFace, with scale = ', num2str(scale)));
%     figure; imshow(result(i).rotationFace); title(strcat('toneFace, with rotation = ', num2str(rotation), ' degrees'));
    perc = uint8(100*(i/length(result)));
    waitbar(i/length(result), h, sprintf('Done %d%%', perc));
end
toc
close(h)











