function [leftEyeCenterX, leftEyeCenterY, ...
          rightEyeCenterX, rightEyeCenterY, ...
          rgbImage, face, mouthMap] = rotateFace(leftEyeCenterX, leftEyeCenterY, ...
                                                 rightEyeCenterX, rightEyeCenterY, ...
                                                 rgbImage, face, mouthMap)

    referenceDirection = [0 1];
    
    eyeDirection = [leftEyeCenterX-rightEyeCenterX, leftEyeCenterY-rightEyeCenterY];
    eyeDirection = eyeDirection / norm(eyeDirection);
    
    angle = acos(dot(referenceDirection, eyeDirection)) * 180 / pi;
    angle = 90 - angle;
    
    marker=zeros(size(face(:,:,1)));
    originalMarker = marker;
    
    rgbImage = imrotate(rgbImage, -angle);
    face = imrotate(face, -angle);
    mouthMap = imrotate(mouthMap, -angle);
    
    marker(round(leftEyeCenterY), round(leftEyeCenterX)) = 1;
    marker_rot = imrotate(marker, -angle, 'bicubic');
    [leftEyeCenterY, leftEyeCenterX]=find(marker_rot > 0);
    leftEyeCenterY = leftEyeCenterY(1);
    leftEyeCenterX = leftEyeCenterX(1);
    
    marker = originalMarker;
    marker(round(rightEyeCenterY), round(rightEyeCenterX)) = 1;
    marker_rot = imrotate(marker, -angle, 'bicubic');
    [rightEyeCenterY, rightEyeCenterX]=find(marker_rot > 0);
    rightEyeCenterY = rightEyeCenterY(1);
    rightEyeCenterX = rightEyeCenterX(1);

end