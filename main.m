function main()
  addpath src

  [access_images, number_of_access_images] = readAllFromDir('access', 'img/access/', '*.jpg');
  [no_access_images, number_of_no_access_images] = readAllFromDir('no_access', 'img/no_access/', '*.jpg');
  [hard_images, number_of_hard_images] = readAllFromDir('hard', 'img/hard/', '*.jpg');


  rmpath src
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read all images from a given directory.
%
% name: Name of output file.
%
% directory: Directory containgin the files.
%
% extension: Extension of the files to read.
%
% images: Matrix containing the images.
%
% number_of_images: Number of images.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [images, number_of_images] = readAllFromDir(name, directory, extension)
  file_names = dir(strcat(directory, extension));
  number_of_images = length(file_names);

  mat_name = strcat(name, '.mat');

  if exist(mat_name)
    images = load(mat_name);
  else
    for k = 1:number_of_images
      current_filename = file_names(k).name;
      current_image = imread(strcat(directory, current_filename));
      images{k} = current_image;
    end

    save(mat_name, 'images');
  end
end
