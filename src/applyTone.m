function output = applyTone(img, lightPercentage)

imgHSV = rgb2hsv(img);
    
if nargin > 1
    V = imgHSV(:,:,3)*lightPercentage;
    
    if lightPercentage > 1
        bin = V > 1;
        V(bin) = 1;
    end
      
    imgHSV(:,:,3) = V;
else
    randNumber = randi([0.8 1]);
    imgHSV(:,:,3) = imgHSV(:,:,3)*randNumber;
end

imgTone = hsv2rgb(imgHSV);
output = uint8(255*imgTone);