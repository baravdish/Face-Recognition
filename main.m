clear
clc
close all

addpath src

[access_images, number_of_access_images] = readAllFromDir('access', 'img/access/', '*.jpg');
[no_access_images, number_of_no_access_images] = readAllFromDir('no_access', 'img/no_access/', '*.jpg');
[hard_images, number_of_hard_images] = readAllFromDir('hard', 'img/hard/', '*.jpg');

ASA = tnm034(access_images{7});

% ASA = tnm034(hard_images{3});

BSB = tnm034(hard_images{3});
% BSB = tnm034(hard_images{4});

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
A = A(1:minLength);
B = B(1:minLength);

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
% figure % opens new figure window
% plot(x, y1Filt, x, y2Filt);
% title('Filtered Comparison');

% hsv = norm(sqrt(abs(AHSV(1,1)-BHSV(1,1))^2 + abs(AHSV(1,2)-BHSV(1,2))^2 + abs(AHSV(1,3)-BHSV(1,3))^2));
% corre(1,2) > 0.9 &&
% Argb = (Argb - min(Argb))./(max(Argb) - min(Argb));
% Brgb = (Brgb - min(Brgb))./(max(Brgb) - min(Brgb));
Argb
Brgb
 iris = min([abs(Argb(1)-Brgb(1)), ...
             abs(1+Argb(1)-Brgb(1)), ...
              abs(Argb(1)-Brgb(1)-1)])
% iris = sum(abs(Argb(1) - Brgb(1)))
hit = iris < 0.2 && hsv > 0.82 && rmse > 0.8 && like > 0.78 && euc > 0.94 && meanPctDiff > 0.85