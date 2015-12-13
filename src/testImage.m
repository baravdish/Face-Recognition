function id = testImage(img, database, size_of_database)
  img = rgb2gray(img);
  img = imresize(img, [200, 200]);
  
  histogram = describe(img, 21);
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

  threshold = 0.001140;
  if(value < threshold)
    id = index;
  else
    id = 0;
  end

end