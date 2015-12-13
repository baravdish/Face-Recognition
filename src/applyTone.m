function output = applyTone(img, lightPercentage)
    
if nargin > 1
    imgLight = img*lightPercentage;
    if lightPercentage > 1
        index = imgLight > 255;
        imgLight(index) = 255;
    end    
else
    randNumber = randi([0.7 1]);
    imgLight = img*randNumber;
end

output = imgLight;