function [] = stats(build_database, write_output, break_at_failure)
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
  lowest_no_access_comparison = realmax;
  height = 200; % Requires that the database is rebuilt!
  width = 200; % Requires that the database is rebuilt!
  threshold = 0.001240; % Requires that the database is rebuilt!
  kernel_size = 21; % Requires that the database is rebuilt!
  decorr = 0; % Requires that the database is rebuilt!
  freqestim = 1; % Requires that the database is rebuilt!

  if build_database == 1
    createDatabase(height, width, kernel_size, decorr, freqestim);
  end

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
      % imshow(image{1})
      has_succeded = false;
      tic
      [id, id_false, min_value] = tnm034(image{1}, height, width, threshold, kernel_size, decorr, freqestim);
      time(n, k) = toc;

      if n == number_of_people
        if lowest_no_access_comparison > min_value
          lowest_no_access_comparison = min_value;
        end
        if id == 0
          success_counter = success_counter + 1;
          has_succeded = true;
        else
          fail_counter = fail_counter + 1;
          false_positives = false_positives + 1;
          has_succeded = false;
          if break_at_failure == 1
            imshow(image{1});
            clc
            disp('False Positive!');
            disp(sprintf('Should be %i. \nWas %i. \nBelived to be %i', n, id, id_false));
            disp(sprintf('Comparision value was %f with threshold %f.', min_value, threshold));
            pause
          end
        end
      else
        if id == n
          success_counter = success_counter + 1;
          has_succeded = true;
        else
          fail_counter = fail_counter + 1;
          has_succeded = false;


          if id == 0
            false_negatives = false_negatives + 1;
          else
            failed_recognitions = failed_recognitions + 1;
          end

          if break_at_failure == 1
            imshow(image{1});
            clc
            if id == 0
              disp('False Negative!');
            else
              disp('Failed Recognition!');
            end
            disp(sprintf('Should be %i. \nWas %i. \nBelived to be %i', n, id, id_false));
            disp(sprintf('Comparision value was %f with threshold %f.', min_value, threshold));
            pause
          end
        end
      end
      images_computed = images_computed + 1;
      if has_succeded == true
        progress(images_computed) = 'T';
      else
        progress(images_computed) = 'F';
      end
      clc
      disp(progress);
    end
  end

  clc
  disp(sprintf('Number of images: %i', total_number_of_images));
  disp(sprintf('Number of successes: %i', success_counter));
  disp(sprintf('Number of failures: %i', fail_counter));
  disp(sprintf('Number of false positives: %i', false_positives));
  disp(sprintf('Number of false negatives: %i', false_negatives));
  disp(sprintf('Number of failed recognitions: %i', failed_recognitions));
  disp(sprintf('Lowest no access comparision: %f', lowest_no_access_comparison));

  fail_rate = fail_counter / (success_counter + fail_counter);
  success_rate = success_counter / (success_counter + fail_counter);

  disp(sprintf('Success rate: %i%%', round(success_rate * 100)));
  disp(sprintf('Fail rate: %i%%', round(fail_rate * 100)));

  average_time = 0;
  total_time = 0;
  shortest_time = realmax;
  longest_time = 0;
  for n = 1 : number_of_people
    for k = 1 : number_of_images
      average_time = average_time + time(n, k);
      if time(n, k) > longest_time
        longest_time = time(n, k);
      end

      if time(n, k) < shortest_time
        shortest_time = time(n, k);
      end
    end
  end

  total_time = average_time;
  average_time = average_time / double(total_number_of_images);
  disp(sprintf('Total time: %0.2fs', total_time));
  disp(sprintf('Average time: %0.2fs', average_time));
  disp(sprintf('Longest time: %0.2fs', longest_time));
  disp(sprintf('Shortest time: %0.2fs', shortest_time));

  if write_output == 1
    format_spec = ['==========Parameters==========\n' ...
                  'Height: %i\n' ...
                  'Width: %i\n' ...
                  'Threshold: %f\n' ...
                  'Kernel size: %i\n' ...
                  'Decorrelation: %i\n' ...
                  'Local frequency estimation: %i\n' ...
                  '==========Results==========\n' ...
                  'Number of images: %i\n' ...
                  'Number of successes: %i\n' ...
                  'Number of failures: %i\n' ...
                  'Number of false positives: %i\n' ...
                  'Number of false negatives: %i\n' ...
                  'Number of failed recognitions: %i\n' ...
                  'Lowest no access comparision: %f\n' ...
                  'Success rate: %i%%\n' ...
                  'Fail rate: %i%%\n' ...
                  'Total time: %0.2fs\n' ...
                  'Average time: %0.2fs\n' ...
                  'Longest time: %0.2fs\n' ...
                  'Shortest time: %0.2fs\n' ...
                  '==============================\n\n'];
    fid = fopen('stats.txt', 'a');
    fprintf(fid, format_spec, height, width, threshold, kernel_size, decorr, freqestim, total_number_of_images, success_counter, fail_counter, false_positives, false_negatives, failed_recognitions, lowest_no_access_comparison, round(success_rate * 100), round(fail_rate * 100), total_time, average_time, longest_time, shortest_time);
  end
end
