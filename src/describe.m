%Inspired by the work of Janne Heikkil, Esa Rahtu, and Ville Ojansivu.
%http://www.cse.oulu.fi/CMV/Downloads/LPQMatlab

function descriptor = describe(img, kernel_size)
  img = double(img);
  radius = (kernel_size - 1) / 2;
  x = -radius : radius;

  w0 = ones(size(x));
  w1 = exp(-2 * pi * x * i / kernel_size);
  w2 = real(w1) - i * imag(w1);

  filter_response = conv2(conv2(img, w0.', 'valid'), w1, 'valid');
  frequency_response = zeros(size(filter_response, 1), size(filter_response, 2), 8);
  frequency_response(:, :, 1) = real(filter_response);
  frequency_response(:, :, 2) = imag(filter_response);

  filter_response = conv2(conv2(img, w1.', 'valid'), w0, 'valid');
  frequency_response(:, :, 3) = real(filter_response);
  frequency_response(:, :, 4) = imag(filter_response);

  filter_response = conv2(conv2(img, w1.', 'valid'), w1, 'valid');
  frequency_response(:, :, 5) = real(filter_response);
  frequency_response(:, :, 6) = imag(filter_response);

  filter_response = conv2(conv2(img, w1.', 'valid'), w2, 'valid');
  frequency_response(:, :, 7) = real(filter_response);
  frequency_response(:, :, 8) = imag(filter_response);

  [rows, columns, number_of_components] = size(frequency_response);

  descriptor = zeros(rows, columns);
  for k = 1 : number_of_components
    binary_representation = frequency_response(:, :, k) > 0;
    descriptor = descriptor + binary_representation * (2^(k - 1));
  end

  descriptor = hist(descriptor(:), 0 : 255);
  descriptor = descriptor / sum(descriptor);

end
