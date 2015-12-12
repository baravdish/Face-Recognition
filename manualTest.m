clear
clc

load resultAccess.mat
% load resultNoAcc.mat
% load resultHard.mat
% load resultTest.mat

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





