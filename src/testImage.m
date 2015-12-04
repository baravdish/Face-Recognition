function [id, id_false, min_value, threshold] = testImage(img, database, size_of_database)
  % img = equalize(img);
  img = rgb2gray(img);
  img = imresize(img, [200, 200]);
  % img = histeq(img);
  % img = LogAbout(img);

  histogram = lpq(img, 21, 1, 2);

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

  threshold = 0.001350;
  if(value < threshold)
    id = index;
  else
    id = 0;
  end
  min_value = value;

end
