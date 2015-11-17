function [img, meanColor, modeColor] = getBackgroundColor (img)
%
% function [img, meanColor, modeColor] = getBackgroundColor (img)
%
%    img   -   Either a string representing the filename of an image to open
%              or an image itself.  If the latter, it must be either a
%              3-dimensional matrix representing an RGB image or a 2-dimensional
%              matrix representing a grayscale image.

if ischar(img)
  img = imread(imageFile);
end
img = double(img);

% Handle RGB and Grayscale separately.
if ndims(img)==3
  % There are probably some spiffy ways to consolidate this sprawl
  % so that the R, G, and B channels are not being processed
  % independently, but for the time being, this does work.
  red   = getBG(img(:, :, 1));
  green = getBG(img(:, :, 2));
  blue  = getBG(img(:, :, 3));

  % For each channel, remove the "foreground" regions identified in
  % each of the other channels.
  red(isnan(green)) = NaN;
  red(isnan(blue)) = NaN;

  green(isnan(red)) = NaN;
  green(isnan(blue)) = NaN;

  blue(isnan(red)) = NaN;
  blue(isnan(green)) = NaN;

  % Compute the mean and mode colors.
  meanColor = [ ...
      mean(mean( red(~isnan(red)) )) ...
      mean(mean( green(~isnan(green)) )) ...
      mean(mean( blue(~isnan(blue)) )) ];
  modeColor = [ ...
      mode(mode( red(~isnan(red)) )) ...
      mode(mode( green(~isnan(green)) )) ...
      mode(mode( blue(~isnan(blue)) )) ];

  % Update each the foreground regions of each channel and set them
  % to their mean colors.  This is only necessary for visualization.
  red(isnan(red)) = meanColor(1);
  green(isnan(green)) = meanColor(2);
  blue(isnan(blue)) = meanColor(3);

  img(:, :, 1) = red;
  img(:, :, 2) = green;
  img(:, :, 3) = blue;
else
  img = getBG(img);
  meanColor = mean(mean( img( ~isnan(img) ) ));
  modeColor = mode(mode( img( ~isnan(img) ) ));
  img(isnan(img)) = meanColor;
end

% Convert the image back to integers (optional)
img = uint8(img);

% Display the results before returning
display(meanColor)
display(modeColor)



  function image = getBG (image)
      mask = getAttenuationMask(size(image), min(size(image)) / 2, 0, 1);
%       figure;
%       imshow(mask);
%       pause;

      % Assume that the background is mostly constant, so isolate the high-frequency
      % parts of the image in the frequency domain and then transform it back into the spatial domain
      fftImage = fftshift(fft2(image));
      fftImage = fftImage .* mask;
      invFftImage = abs(ifft2(fftImage));

      % Expand the high-frequency areas of the image and fill in any holes.  This should
      % cover all but the (hopefully) low frequency background areas.
      edgeRegion = imfill(imdilate(invFftImage, strel('disk', 4, 4)), 'holes');

      % Now remove the parts of the image that are covered by edgeRegion
      edgeMean = mean(mean(edgeRegion));
      image(edgeRegion>edgeMean) = NaN;
  end
end