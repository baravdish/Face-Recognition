%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: RGB-Image to correct.
%
% output: Corrected image.
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = whiteBalance(img)
%LIGHTNORM Summary of this function goes here
%   Detailed explanation goes here

R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);

rVec = R(:);
gVec = G(:);
bVec = B(:);

sortedImg(:,1) = sort(rVec, 'descend');
sortedImg(:,2) = sort(gVec, 'descend'); 
sortedImg(:,3) = sort(bVec, 'descend');

nPixels = length(rVec);
percentage = 0.001;
nFind = round(nPixels*percentage);

meanValues = [mean(sortedImg((1:nFind),1)), ...
              mean(sortedImg((1:nFind),2)), ...
              mean(sortedImg((1:nFind),3))];


scaleFactors = 255./meanValues;


rNorm = R*scaleFactors(1);
binR = rNorm > 255;
rNorm(binR) = 255;

gNorm = G*scaleFactors(2);
binG = gNorm > 255;
gNorm(binG) = 255;

bNorm = B*scaleFactors(3);
binB = bNorm > 255;
bNorm(binB) = 255;
output = cat(3, rNorm, gNorm, bNorm);

% figure
% imshow(output,[])

end