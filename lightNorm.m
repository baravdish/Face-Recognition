function [ result ] = lightNorm( img )
%LIGHTNORM Summary of this function goes here
%   Detailed explanation goes here

imgMean = mean(img, 3);

pixVec = imgMean(:);

sortPixVec = sort(pixVec, 'descend');

nTop = length(sortPixVec)*0.05;

topBright = sortPixVec(1:nTop);

R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);

AR = mean(R(:));
AG = mean(G(:));
AB = mean(B(:));

avgGray = (AR + AG + AB)/3;

aR = avgGray/AR;

aG = avgGray/AG;

aB = avgGray/AB;


r = aR*img(:,:,1);
binR = r > 255;

r(binR) = 255;

g = aG*img(:,:,2);
binG = g > 255;

g(binG) = 255;

b = aB*img(:,:,3);
binB = b > 255;


b(binB) = 255;


result = cat(3,r,g,b);


end

