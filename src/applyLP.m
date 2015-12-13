function output = applyLP(img, sigma)

if nargin > 1 
    output = imgaussfilt(img, sigma);
else
    randNumber = randi([0 1.5]);
    output = imgaussfilt(img, randNumber);
end