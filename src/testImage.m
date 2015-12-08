function [id, id_false, min_value] = testImage(img, database, size_of_database, height, width, threshold, kernel_size, decorr, freqestim)
  % img = equalize(img);
  img = rgb2gray(img);
%   img = imgradient(img);
%   img = imadjust(img);
%   img = imsharpen(img);
% img = histeq(img);


  img = imresize(img, [height, width]);
  % img = histeq(img);
  % img = LogAbout(img);

  histogram = lpq(img, kernel_size, decorr, freqestim);

  hist_comp = zeros(size_of_database, 1);
  small_comp = zeros(size_of_database, 1);

  table = zeros(size_of_database, 256);

  for i = 1 : size_of_database
    for j = 1 : 256
      table(i, j) = abs(database(i, j) - histogram(j));
      hist_comp(i) = hist_comp(i) + table(i, j);
    end
    hist_comp(i) = hist_comp(i) / 256;

  [value index] = min(hist_comp);

  id_false = index;

  if(value < threshold)
    id = index;
  else
    id = 0;
  end
  min_value = value;

end
