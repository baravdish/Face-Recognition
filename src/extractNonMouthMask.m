function nonMouthMask = extractNonMouthMask(face, eyeDistance, eyesMeanY, ...
                                            leftEyeCenterX, rightEyeCenterX, ...
                                            offsetX)
   
    nonMouthMask = zeros(size(face(:,:,1)));
    
    measureDistance = eyeDistance * 0.8;
    
    nonMouthMask(1: eyesMeanY+measureDistance, :) = 1;
    nonMouthMask(1: end, 1:max(leftEyeCenterX-offsetX, 1)) = 1;
    nonMouthMask(1: end, min(rightEyeCenterX+offsetX, end):end) = 1;
    
end