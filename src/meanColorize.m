function binaryImage = meanColorize(rgbImage, mask, tolerance)
    rows = size(rgbImage, 1);
    columns = size(rgbImage, 2);
    
    if nargin < 2
        level = graythresh(rgbImage);
        bw = im2bw(rgbImage, level);
        mask = bw;
    end
    
    if nargin < 3 
        tolerance = 27;
    end

    % Convert image from RGB colorspace to lab color space.
    cform = makecform('srgb2lab');
%     lab_Image = applycform(im2double(rgbImage) ,cform);
%     lab_Image = rgb2lab(rgbImage);
%     lab_Image = rgb2hsv(rgbImage);
    lab_Image = rgb2ycbcr(im2double(rgbImage));

    % Extract out the color bands from the original image
    % into 3 separate 2D arrays, one for each color component.
    LChannel = lab_Image(:, :, 1); 
    aChannel = lab_Image(:, :, 2); 
    bChannel = lab_Image(:, :, 3); 
    
    % Compute mean Lab values
    LVector = LChannel(mask); 
    LMean = mean(LVector(:));
    aVector = aChannel(mask);
    aMean = mean(aVector(:));
    bVector = bChannel(mask);
    bMean = mean(bVector(:));

    % Make uniform images of only that one single LAB color.
    LStandard = LMean * ones(rows, columns);
    aStandard = aMean * ones(rows, columns);
    bStandard = bMean * ones(rows, columns);

    % Create the delta images: delta L, delta A, and delta B.
    deltaL = LChannel - LStandard;
    deltaa = aChannel - aStandard;
    deltab = bChannel - bStandard;
    
    
    % Create the Delta E image.
    % This is an image that represents the color difference.
    % Delta E is the square root of the sum of the squares of the delta images.
    deltaE = sqrt(deltaL.^2 + deltaa.^2 + deltab.^2);

    % Mask it to get the Delta E in the mask region only.
    maskedDeltaE = deltaE .* mask;
    % Get the mean delta E in the mask region
    % Note: deltaE(mask) is a 1D vector of ONLY the pixel values within the masked area.
    meanMaskedDeltaE = mean(deltaE(mask));
    % Get the standard deviation of the delta E in the mask region
    stDevMaskedDeltaE = std(deltaE(mask));
    
%     tolerance = meanMaskedDeltaE + 3 * stDevMaskedDeltaE;
%     tolerance = 27
    
    binaryImage = (deltaE >= meanMaskedDeltaE - tolerance) & ...
                  (deltaE <= meanMaskedDeltaE + tolerance);
%     binaryImage = (deltaE >= meanMaskedDeltaE - 10) & ...
%                   (deltaE <= meanMaskedDeltaE + 10);
%     
    
    
%     matchingColors = bsxfun(@times, rgbImage, cast(binaryImage, class(rgbImage)));
    
%     fontSize = 12;
%     
%     figure
%     imshow(rgbImage)
%     caption = sprintf('Input');
% 	title(caption, 'FontSize', fontSize);
%     
%     figure
%     imshow(mask)
%     caption = sprintf('Mask');
% 	title(caption, 'FontSize', fontSize);
%     
%     figure
%     imshow(matchingColors);
% 	caption = sprintf('Matching Colors (Delta E <= %.1f)', tolerance);
% 	title(caption, 'FontSize', fontSize);
% 
%     nonMatchingColors = bsxfun(@times, rgbImage, cast(~binaryImage, class(rgbImage)));
%     
%     figure
%     imshow(nonMatchingColors);
% 	caption = sprintf('Non-Matching Colors (Delta E > %.1f)', tolerance);
% 	title(caption, 'FontSize', fontSize);
    
end
