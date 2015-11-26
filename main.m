
clc
close all

addpath src

[access_images, number_of_access_images] = readAllFromDir('access', 'img/access/', '*.jpg');
[no_access_images, er_of_no_access_images] = readAllFromDir('no_access', 'img/no_access/', '*.jpg');
[hard_images, number_of_hard_images] = readAllFromDir('hard', 'img/hard/', '*.jpg');

result = tnm034(no_access_images{2});
result = tnm034(no_access_images{4});
result = tnm034(access_images{8});
result = tnm034(hard_images{37});
% result = tnm034(access_images{8});

result = tnm034(no_access_images{1});


result = tnm034(no_access_images{1});
result = tnm034(access_images{8});



result = tnm034(no_access_images{4});
% pause;

result = tnm034(hard_images{5});
% pause
result = tnm034(access_images{6});




result = tnm034(hard_images{36});
result = tnm034(hard_images{33});
result = tnm034(hard_images{15});
result = tnm034(no_access_images{1});
result = tnm034(no_access_images{2});
result = tnm034(no_access_images{3});
result = tnm034(no_access_images{4});
result = tnm034(hard_images{37});
result = tnm034(hard_images{17});
result = tnm034(hard_images{35});

result = tnm034(no_access_images{2});
result = tnm034(hard_images{28});
result = tnm034(no_access_images{1});
% result = tnm034(hard_images{2});
result = tnm034(no_access_images{4});
% result = tnm034(no_access_images{1});
% result = tnm034(hard_images{34});
result = tnm034(hard_images{4});
% result = tnm034(no_access_images{1});
% result = tnm034(hard_images{27});
% result = tnm034(access_images{8});
% result = tnm034(access_images{6});


% result = tnm034(hard_images{34});
% 
% 
% result = tnm034(no_access_images{1});
% result = tnm034(no_access_images{2});


        
% result = tnm034(no_access_images{4});

result = tnm034(no_access_images{1});
result = tnm034(hard_images{28});
result = tnm034(no_access_images{4});
result = tnm034(hard_images{34});


% % result = tnm034(hard_images{15});
% % result = tnm034(no_access_images{2});
% % 
% % result = tnm034(access_images{2});
% % result = tnm034(hard_images{14});
% % result = tnm034(hard_images{11});
% % 
% % result = tnm034(hard_images{33});
% % result = tnm034(hard_images{34});
% result = tnm034(hard_images{37});
% result = tnm034(hard_images{7});
% 
% result = tnm034(hard_images{2});
result = tnm034(no_access_images{1});
result = tnm034(no_access_images{2});
result = tnm034(no_access_images{3});
