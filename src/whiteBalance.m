%%%%%%%%%%%%%%%%%%%%%%%%%
% input: RGB-Image to correct.

% output: Corrected image.
%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = whiteBalance(img)
%LIGHTNORM Summary of this function goes here
%   Detailed explanation goes here
    
    originalImg = img;
    
 
      
%     end
    
%     X = rgb2gray(img);
%     figure
%     imshow(X)
%     title('Original','FontWeight','bold')
%     prev = otsu(X,2);
%     for n = 3:10
%       IDX = otsu(X,n);
%       IDX = IDX - prev;
% %       IDX = imadjust(IDX);
%       figure;
% %       level = graythresh(IDX);
% %       im2bw(IDX, level)
%       imagesc(IDX);
%       title(['n = ' int2str(n)]);
%       pause;
%       prev = IDX;
%     end
%     pause;
    
   
%     figure
%     level = graythresh(img);
%     bw = im2bw(img, level);
% %     bw = bwareaopen(bw, 50);
% %     bw = ~bw;
% %     bw = imfill(bw, 'holes');
%     imshow(bw)
% 
%     pause;
%     figure;
%     overSaturated = rgb2gray(img);
%     %     imshow(overSaturated); title('Over saturation'); pause;
%     %     imshow(imadjust(overSaturated)); title('imadjust(overSaturated)'); pause;
% %     imshow(imfilter(imadjust(overSaturated), fspecial('laplacian'))); title('filter imadjust(overSaturated)'); pause;
%     avg = ones(5,5) / 25;
%     overSaturated = imfilter(overSaturated, avg);
%     overSaturated = imfilter(overSaturated, fspecial('gaussian'));
%     overSaturated = imfilter(overSaturated, fspecial('laplacian'));
%     %     imshow(overSaturated); title('filter Over saturation'); pause;
%     overSaturated = imdilate(overSaturated, strel('disk', 2));
% %     imshow(overSaturated); title('dilate Over saturation'); pause;
%     overSaturated = ~overSaturated;
%     imshow(overSaturated); title('filter Over saturation'); pause;
%     
%     imshow(img); title('before Over saturation'); pause;
%     overSaturated = ~overSaturated;
%     overSaturatedRep = repmat(overSaturated, [1,1,3]);
%     img = img.*uint8(overSaturatedRep);
%     imshow(img); title('after Over saturation'); pause;
%     
%     [pixelCounts grayLevels] = imhist(rgb2gray(originalImg));
%     cdf = cumsum(pixelCounts) / sum(pixelCounts);
%     subplot(1,4, 1);
%     imshow(originalImg);
%     subplot(1,4, 2);
%     bar(pixelCounts);
%     subplot(1,4, 3);
%     plot(cdf);
%     subplot(1,4, 4);
%     thresholdIndex = find(cdf < 0.99, 1, 'last');
%     thresholdValue = grayLevels(thresholdIndex)
%     thresholdedImage = originalImg;
%     thresholdedImage(originalImg(:) < thresholdValue) = 0;
%     imshow(thresholdedImage);
%     set(gcf, 'Position', get(0, 'ScreenSize')); % Maximize figure.
%     pause;
    
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


R = originalImg(:,:,1);
G = originalImg(:,:,2);
B = originalImg(:,:,3);

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

% figure; imshow(originalImg); title('original'); pause;
% figure; imshow(output); title('WhiteBalanced'); pause;

end


% 
% function out = whiteBalance(in)
%     method = 'gray'
% % WHITEBALANCE - white balance an image using gray world or normalization
% % out = whiteBalance(in, method)
% 
% %% -------------------------------------------------------------------------
% % Matlab code and data to reproduce results from the paper                  
% % "Joint Demosaicing and Super-Resolution Imaging from a Set of             
% % Unregistered Aliased Images"                                              
% % Patrick Vandewalle, Karim Krichane, David Alleysson and Sabine S?sstrunk  
% % available at http://lcavwww.epfl.ch/reproducible_research/VandewalleKAS07/
% %                                                                           
% % Copyright (C) 2007 Laboratory of Audiovisual Communications (LCAV),       
% % Ecole Polytechnique Federale de Lausanne (EPFL),                          
% % CH-1015 Lausanne, Switzerland.                                            
% %                                                                           
% % This program is free software; you can redistribute it and/or modify it   
% % under the terms of the GNU General Public License as published by the     
% % Free Software Foundation; either version 2 of the License, or (at your    
% % option) any later version. This software is distributed in the hope that  
% % it will be useful, but without any warranty; without even the implied     
% % warranty of merchantability or fitness for a particular purpose.          
% % See the GNU General Public License for more details                       
% % (enclosed in the file GPL).                                               
% %                                                                           
% % Latest modifications: June 7, 2007.                                       
% 
%     if nargin == 1
%         method = 'gray';
%     end
%     
%     switch(method)
%         case 'norm' % Channel normalization
%             out(:,:,1) = in(:,:,1) * (1/max(max(in(10:end-10,10:end-10,1))));
%             out(:,:,2) = in(:,:,2) * (1/max(max(in(10:end-10,10:end-10,2))));
%             out(:,:,3) = in(:,:,3) * (1/max(max(in(10:end-10,10:end-10,3))));
%         case 'gray' % Gray world
%             mG = mean(mean(in(10:end-10,10:end-10,2)));
%             mR = mean(mean(in(10:end-10,10:end-10,1)));
%             mB = mean(mean(in(10:end-10,10:end-10,3)));
%             out(:,:,1) = in(:,:,1) * (mG/mR);
%             out(:,:,3) = in(:,:,3) * (mG/mB);
%             out(:,:,2) = in(:,:,2);
%             out = out / max(max(max(out(10:end-10,10:end-10, :))));
%         otherwise
%             error('unknown method')
%     end
% end