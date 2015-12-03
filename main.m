clc
close all

addpath src

[access_images, number_of_access_images] = readAllFromDir('access', 'img/access/', '*.jpg');
[no_access_images, er_of_no_access_images] = readAllFromDir('no_access', 'img/no_access/', '*.jpg');
[hard_images, number_of_hard_images] = readAllFromDir('hard', 'img/hard/', '*.jpg');

im = hard_images{15};
imshow(im)
result = tnm034(im)
return;
%%
for i = 1:length(access_images)
   
    result = tnm034(access_images{i})

end

%%
for i = 1:length(hard_images)
    
    result = tnm034(hard_images{i})
    
end

%% 
for i = 1:length(no_access_images)
    
    result = tnm034(hard_images{i})
    
end
