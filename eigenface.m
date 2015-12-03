clear
clc
close all

addpath src

[access_images, number_of_access_images] = readAllFromDir('access', 'img/cut/', '*.png');

%A = [200 000, antal inlästa bilder]
A = zeros(500*400, number_of_access_images);
average_face = zeros(500*400, 1);
image_vector = zeros(500*400, number_of_access_images);
for k = 1 : number_of_access_images
  image_gray = rgb2gray(access_images{k});
  resized_image = imresize(image_gray, [500 400]);
  image_vector(:, k) = double(reshape(resized_image, [500*400 1])) / 255;
  average_face  = average_face + image_vector(:, k);
end
%normalisera average_face
average_face = average_face / number_of_access_images;

%colormap gray
%imagesc(reshape(average_face, [500, 400]));

for k = 1 : number_of_access_images
  phi = image_vector(:,k) - average_face;
  A(:, k) = phi;
end

L = A'*A;

%V är en matris med egenvektorer och D är en diagonal matris med egenvärdena
[V,D] = eig(L);

eigenfaces = zeros(500*400, number_of_access_images);

figure
colormap gray
for l = 1 : number_of_access_images
  for k = 1 : number_of_access_images
    %ekv. (6) artikel Eigenfaces_Turk sida 75
    eigenfaces(:, l) = eigenfaces(:, l) + V(l,k) * A(:,k);
  end
  subplot(4,4,l)
  imagesc(reshape(eigenfaces(:, l), [500, 400]));
end

known_omegas = zeros(number_of_access_images, number_of_access_images);
for l = 1 : number_of_access_images
  for k = 1 : number_of_access_images
    %ekv. (7) artikel Eigenfaces_Turk sida 75
    known_omegas(k, l) = eigenfaces(:, k)' * (image_vector(:, l) - mean(image_vector(:, l)));
  end
end

% Verify
image_to_check = imread('img/cut_test/4.png');
image_to_check_gray = rgb2gray(image_to_check);
figure
colormap gray
imagesc(image_to_check_gray);
image_to_check_resized = imresize(image_to_check_gray, [500 400]);
image_to_check_vector = double(reshape(image_to_check_resized, [500*400 1])) / 255;
unknown_omega = zeros(number_of_access_images, 1);

figure
colormap gray
for k = 1 : number_of_access_images
  unknown_omega(k) = eigenfaces(:, k)' * (image_to_check_vector - mean(image_to_check_vector));
  subplot(4,4,k)
  imagesc(reshape(image_vector(:, k), [500, 400]));
end


for k = 1 : number_of_access_images
  %ekv. (8) artikel Eigenfaces_Turk sida 76
  distances(k) = norm(unknown_omega - known_omegas(:, k));
end

%printar ut distansen, lågt värde = ansiktet tillhör databasen
distances
