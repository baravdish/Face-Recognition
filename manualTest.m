clear
clc

% load resultAccess.mat
% load resultAccess2.mat
% load resultAccess3.mat
% load resultAccessLight.mat
% load resultAccessLight2.mat

% load resultNoAcc.mat
% load resultNoAcc2.mat
% load resultNoAcc3.mat
% load resultNoAccLight.mat
% load resultNoAccLight2.mat

% load resultHard.mat
% load resultHard2.mat
% load resultHard3.mat
% load resultHardLight.mat
% load resultHardLight2.mat

% load resultTest.mat
% load resultTest2.mat
% load resultTest3.mat
% load resultTestLight.mat
% load resultTestLight2.mat
% load resultTestRot.mat
load resultTestRot2.mat
result(1)
%% Evaluation
% Blur
close all

for i = 1:length(result)

    figure; 
    imshow(result(i).blurFace); 
    title(strcat('blurFace, with sigma = ', num2str(result(i).sigma)));
end

%% Tone
close all

for i = 1:length(result)
   
    figure; imshow(result(i).toneFace);
    title(strcat('toneFace, with tonevalue = ', num2str(result(i).light)));
end

%% Scale
close all

for i = 1:length(result)
   
    figure;
    imshow(result(i).scaleFace); 
    title(strcat('scaleFace, with scale = ', num2str(result(i).scale)));
end

%% Rotation
close all

for i = 1:length(result)
    
    figure; 
    imshow(result(i).rotationFace); 
    title(strcat('rotationFace, with rotation = ', num2str(result(i).rotation), ' degrees'));

end





