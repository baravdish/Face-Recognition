function [] = stats()
  clc
  close all
  warning off
  addpath src

  image_collection = struct;
  number_of_people = 17;
  success_counter = 0;
  fail_counter = 0;
  total_number_of_images = 0;
  false_positives = 0;
  false_negatives = 0;
  failed_recognitions = 0;
  fatal_errors = 0;

  for n = 1 : number_of_people
    folder_name = sprintf('img/tests/%d/', n);
    mat_name = sprintf('p%d', n);
    [images, number_of_images] = readAllFromDir(mat_name, folder_name, '*.jpg');
    image_collection.images{n} = images;
    image_collection.number_of_images{n} = number_of_images;
    total_number_of_images = total_number_of_images + number_of_images;
  end
 
    
  progress = '';
  for n = 1 : total_number_of_images
    progress = strcat(progress, '-');
  end
  
 
  images_computed = 0;
  for n = 1 : number_of_people
    images = image_collection.images{n};
    number_of_images = image_collection.number_of_images{n};

    for k = 1 : number_of_images
      image = images(k);
      id = tnm034(image{1});
      if n == number_of_people
        if id == 0
          success_counter = success_counter + 1;
        else
          fail_counter = fail_counter + 1;
        end
      else
        if id == n
          success_counter = success_counter + 1;
        else
          fail_counter = fail_counter + 1;
        end
      end
      images_computed = images_computed + 1;
      clc
      disp(progress);
    end
  end

  success_rate = success_counter / (images_computed);
  
  disp(sprintf('Success rate: %i%%', round(success_rate * 100)));
end
