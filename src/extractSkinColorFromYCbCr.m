function [skinMask, Y, Cb, Cr] = extractSkinColorFromYCbCr(rgbImage)
    
    imgYCbCr = rgb2ycbcr(rgbImage);
    
    Y = imgYCbCr(:,:,1);
    Cb = imgYCbCr(:,:,2);
    Cr = imgYCbCr(:,:,3);
    
    % Values from:
    % https://web.stanford.edu/class/ee368/Project_03/Project/reports/ee368group15.pdf
    
    skinMask = and(and(Cb > 100, Cb < 145), ... 
                   and(Cr > 132, Cr < 165));
    
%     figure; imshow(Y); title('Y'); pause;
%     figure; imshow(Cb); title('Cb'); pause;
%     figure; imshow(Cr); title('Cr'); pause;
    
end