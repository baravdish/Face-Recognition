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
