clear
clc

% TODO: L�S IN DB. PRE-COMPUTED
% load mat

nImages = 4;
path = 'DB0\db0_';
fileformat = '.jpg';

for i = 1:nImages
    
  img = imread(strcat(path, num2str(i), fileformat));
  figure;
  imshow(img);
  
  % TODO: Pipeline f�r face recognition
  
  
  
end

pause
close all