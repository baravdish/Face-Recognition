
clc
close all

addpath src

[access_images, number_of_access_images] = readAllFromDir('access', 'img/access/', '*.jpg');
[no_access_images, er_of_no_access_images] = readAllFromDir('no_access', 'img/no_access/', '*.jpg');
[hard_images, number_of_hard_images] = readAllFromDir('hard', 'img/hard/', '*.jpg');

for i = 1:length(access_images)

out = tnm034(access_images{i});
filename = (['accessImg_',num2str(i),'.png']);
imwrite(out, filename);
end

%%
for i = 1:length(hard_images)

out = tnm034(hard_images{i});
filename = (['hardImg_',num2str(i),'.png']);
imwrite(out, filename);
end

%% 
for i = 1:length(no_access_images)
    out = tnm034(no_access_images{i});
    filename = (['noAccessImg',num2str(i),'.png']);
    imwrite(out, filename);
end