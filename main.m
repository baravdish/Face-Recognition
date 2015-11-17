clear
clc
close all

addpath src
% 
% filename = 'img/access/db1_04.jpg';
% outFile = 'img/access/db1_04_out.jpg';
% option = 'cat';
% catType = 'CAT02';
% maxIter = 1000;
% T = 0.3;
% plot = 1;
% robustAWB(filename,outFile,option,catType,T,maxIter,plot);

[access_images, number_of_access_images] = readAllFromDir('access', 'img/access/', '*.jpg');
[no_access_images, number_of_no_access_images] = readAllFromDir('no_access', 'img/no_access/', '*.jpg');
[hard_images, number_of_hard_images] = readAllFromDir('hard', 'img/hard/', '*.jpg');

    
    
    
    
% I = hard_images{37};
% if(size(I, 3) > 1)
%     final_image = zeros(size(I,1), size(I,2));
%     for i = 1:size(I,1)
%         for j = 1:size(I,2)
%             R = I(i,j,1);
%             G = I(i,j,2);
%             B = I(i,j,3);
%             if(R > 95 && G > 40 && B > 20)
%                 v = [R,G,B];
%                 if((max(v) - min(v)) > 15)
%                     if(abs(R-G) > 15 && R > G && R > B)
%                         %it is a skin
%                         final_image(i,j) = 1;
%                     end
%                 end
%             end
%         end
%     end
% end %%% added
% 
% final_image = ExtractNLargestBlobs(final_image, 1);
% figure, imshow(final_image);
% disp('before');
% BW=final_image;
% % else
% BW=im2bw(I);
% % figure, imshow(BW);
% % disp('convert');
% % end

% %%%%%%%%%%%Connected Components-Detection of face%%%%%%%%%%%%%
% L = bwlabel(BW,4); %BWLABEL Label==>connected components in 2-D binary image
% BB = regionprops(L, 'BoundingBox');%Measure properties of image regions
% BB1 =struct2cell(BB);%struct2cell==>Convert structure array to cell array
% BB2 = cell2mat(BB1);%cell2mat==>Convert a cell array into a single matrix.
% [s1 s2]=size(BB2);%sizing BB2 in matrix format
% mx=0;
% for k=3:4:s2-1
%     p=BB2(1,k)*BB2(1,k+1);
%     if p>mx && (BB2(1,k)/BB2(1,k+1))<1.8
%         mx=p;
%         j=k;
%     end
% end
% 
% return;


% 
% original = imadjust(rgb2gray(access_images{1}));
% 
% scale = 1.3;
% J = imresize(original, scale);
% theta = 45;
% distorted = imrotate(J,theta);
% % distorted = imadjust(rgb2gray(hard_images{33}));
% 
% 
% 
% return
% 
% figure
% imshow(distorted)
% 
% ptsOriginal  = detectSURFFeatures(original);
% ptsDistorted = detectSURFFeatures(distorted);
% 
% [featuresOriginal,validPtsOriginal]  = extractFeatures(original,ptsOriginal);
% [featuresDistorted,validPtsDistorted] = extractFeatures(distorted,ptsDistorted);
% 
% indexPairs = matchFeatures(featuresOriginal,featuresDistorted);
% 
% matchedOriginal  = validPtsOriginal(indexPairs(:,1));
% matchedDistorted = validPtsDistorted(indexPairs(:,2));
% 
% figure
% showMatchedFeatures(original,distorted,matchedOriginal,matchedDistorted)
% title('Candidate matched points (including outliers)')
% 
% [tform, inlierDistorted,inlierOriginal] = estimateGeometricTransform(matchedDistorted,matchedOriginal,'similarity');
% 
% figure
% showMatchedFeatures(original,distorted,inlierOriginal,inlierDistorted)
% title('Matching points (inliers only)')
% legend('ptsOriginal','ptsDistorted')
% 
% outputView = imref2d(size(original));
% recovered  = imwarp(distorted,tform,'OutputView',outputView);
% 
% figure
% imshowpair(original,recovered,'montage')
% % J = hard_images{10};
% % 
% % corners = detectHarrisFeatures(rgb2gray(I));
% % 
% % figure
% % imshow(I); hold on;
% % plot(corners.selectStrongest(50));
% % 
% % corners = detectHarrisFeatures(rgb2gray(J));
% % figure
% % imshow(J); hold on;
% % plot(corners.selectStrongest(50));
% 
% 
% return 
% robustAWB(hard_images{7},hard_images{1});

% figure 
% imshow(hard_images{3});
% 
% figure 
% imshow(whiteBalance(hard_images{3}));
% 
% img = hard_images{3};
% R = img(:,:,1);
% G = img(:,:,2);
% B = img(:,:,3);
% 
% x = double(R);  
% m = mean(x(:));  
% s = std(x(:));
% R = (x - m) / s;  
% 
% x = double(G);  
% m = mean(x(:));  
% s = std(x(:));
% G = (x - m) / s;   
% 
% x = double(B);  
% m = mean(x(:));  
% s = std(x(:));
% B = (x - m) / s;   
% 
% im = img;
% im(:,:,1) = R; im(:,:,2) = G; im(:,:,3) = B;
% 
% figure 
% imshow(im*255);
% 
% R = double(img(:,:,1));
% G = double(img(:,:,2));
% B = double(img(:,:,3));
% 
% NormalizedRed = R(:,:)./sqrt(R(:,:).^2+G(:,:).^2+B(:,:).^2);
% NormalizedGreen = G(:,:)./sqrt(R(:,:).^2+G(:,:).^2+B(:,:).^2);
% NormalizedBlue = B(:,:)./sqrt(R(:,:).^2+G(:,:).^2+B(:,:).^2);
% 
% norm(:,:,1) = NormalizedRed(:,:);
% norm(:,:,2) = NormalizedGreen(:,:);
% norm(:,:,3) = NormalizedBlue(:,:);
% figure 
% imshow(norm)
% 
% figure 
% imshow(color_normalization(img));
% 
% return 

% vec = [1 2 3 4 5 6 7 8 9 10];
% interpolated = imresize(vec,[1 5], 'bicubic')
% % 'nearest'	
% % 'bilinear'	
% % 'bicubic'
% return

rotated = imread('img/rotated.jpg');
% ASA = tnm034(hard_images{37});
% ASA = tnm034(no_access_images{1});

% [img, meanColor, modeColor] = getBackgroundColor(no_access_images{2});
% figure;
% imshow(img);
% pause;

% ASA = tnm034(no_access_images{2});
% 
ASA = tnm034(hard_images{5});
% 

% % ASA = tnm034(hard_images{5});
% 
% % ASA = tnm034(hard_images{1});
% ASA = tnm034(hard_images{2});
% ASA = tnm034(hard_images{3});
% ASA = tnm034(hard_images{4});
% ASA = tnm034(hard_images{5});
% ASA = tnm034(hard_images{6});
% ASA = tnm034(hard_images{7});
% ASA = tnm034(hard_images{8});
% ASA = tnm034(hard_images{9});
% ASA = tnm034(hard_images{10});
% ASA = tnm034(hard_images{11});
% ASA = tnm034(hard_images{12});
% ASA = tnm034(hard_images{13});
% ASA = tnm034(hard_images{14});
% ASA = tnm034(hard_images{15});
% ASA = tnm034(hard_images{16});
% ASA = tnm034(hard_images{17});
% ASA = tnm034(hard_images{18});
% ASA = tnm034(hard_images{19});
% ASA = tnm034(hard_images{20});
% 
% ASA = tnm034(hard_images{2});
% BSB = tnm034(hard_images{7});

% BSB = tnm034(hard_images{9});
% ASA = tnm034(hard_images{1});

% ASA = tnm034(no_access_images{2});
BSB = tnm034(hard_images{1});






A = ASA(1:length(ASA)-7);
B = BSB(1:length(BSB)-7);
SA = ASA(length(ASA)-6);
SB = BSB(length(BSB)-6);
AHSV = ASA(length(ASA)-5:length(ASA)-2)
BHSV = BSB(length(BSB)-5:length(BSB)-2)
Argb = ASA(length(ASA)-2:end);
Brgb = BSB(length(BSB)-2:end);

size(A);
size(B);
minLength = min(SA, SB);
A = A(1:SA);
B = B(1:SB);
maxLength = max(SA, SB);
A = imresize(A,[1 maxLength], 'bicubic');
B = imresize(B,[1 maxLength], 'bicubic');

% A = A(1:minLength);
% B = B(1:minLength);

x = 1:length(A);
y1 = A;
y2 = B;
windowSize = 9;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
y1Filt = filter(b,a,y1);
y2Filt = filter(b,a,y2);

A = y1Filt;
B = y2Filt;


As=(A-min(min(A)))/(max(max(A))-min(min(A)));
Bs=(B-min(min(B)))/(max(max(B))-min(min(B)));

As2 = (A - min(A))./(max(A) - min(A));
Bs2 = (B - min(B))./(max(B) - min(B));
% diff(f1(1:minLength), f2(1:minLength))

euc = 1 - dot(A-B, A-B)/sqrt(dot(A,A)*dot(B,B))
rmse = 1 - sqrt(mean((As-Bs).^2))

like = 1;
% like2 = 0;
for n = 2 : length(A)
    aDiff = A(n)-A(n-1);
    bDiff = B(n)-B(n-1);
    
    
    if (aDiff > 0 && bDiff > 0) || ...
       (aDiff < 0 && bDiff < 0) || ...
       (aDiff == 0 && bDiff == 0)
        like = like + 1;
%         like2 = like2 + abs(As2(n)-Bs2(n));
%         like = like + 1*(length(A)-n);

%         like = like + min(A(n)/B(n), B(n)/A(n));
%         like = like + abs(A(n)-B(n));
    end
%     else
%         like = like + min(A(n)/B(n), B(n)/A(n));
%     end
%     like = like + min(A(n)/B(n), B(n)/A(n));
%     like = like + 1-abs(aDiff-bDiff);
end
% like = (like / length(A))
like = (like / length(A))
% like2

deltaSignal = abs(B - A);
percentageDifference = deltaSignal ./ A; % Percent by element.
meanPctDiff = 1 - mean(percentageDifference) % Average percentage over all elements.

a = min([abs(AHSV(1,1)-BHSV(1,1)), abs(1+AHSV(1,1)-BHSV(1,1)), abs(AHSV(1,1)-BHSV(1,1)-1)]);
b = min([abs(AHSV(1,2)-BHSV(1,2)), abs(1+AHSV(1,2)-BHSV(1,2)), abs(AHSV(1,2)-BHSV(1,2)-1)]);
c = min([abs(AHSV(1,3)-BHSV(1,3)), abs(1+AHSV(1,3)-BHSV(1,3)), abs(AHSV(1,3)-BHSV(1,3)-1)]);
% a
% b
% hsv = 1 - (a+b) / 2
hsv = 1 - (a) / 2

% [corre,P] = corrcoef(A, B)
% hd = HausdorffDist(A,B) 

% figure % opens new figure window
% plot(x, y1, x, y2);
% title('Comparison');

% plot(x,y)

% 
figure % opens new figure window
% plot(x, y1Filt, x, y2Filt);
plot(x, A, x, B);
title('Filtered Comparison');

% hsv = norm(sqrt(abs(AHSV(1,1)-BHSV(1,1))^2 + abs(AHSV(1,2)-BHSV(1,2))^2 + abs(AHSV(1,3)-BHSV(1,3))^2));
% corre(1,2) > 0.9 &&
% Argb = (Argb - min(Argb))./(max(Argb) - min(Argb));
% Brgb = (Brgb - min(Brgb))./(max(Brgb) - min(Brgb));
% Argb
% Brgb
 iris = min([abs(Argb(1)-Brgb(1)), ...
             abs(1+Argb(1)-Brgb(1)), ...
              abs(Argb(1)-Brgb(1)-1)])
% iris = sum(abs(Argb(1) - Brgb(1)))
hit = iris < 0.2 && hsv > 0.82 && rmse > 0.8 && like > 0.78 && euc > 0.94 && meanPctDiff > 0.85