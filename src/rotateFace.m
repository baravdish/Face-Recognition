function [leftEyeCenterX, leftEyeCenterY, ...
          rightEyeCenterX, rightEyeCenterY, ...
          rgbImage, face, faceMask] = rotateFace(leftEyeCenterX, leftEyeCenterY, ...
                                       rightEyeCenterX, rightEyeCenterY, ...
                                       rgbImage, face, faceMask)

%     figure; imshow(face); title('face before rotated'); 
%     viscircles([leftEyeCenterX, leftEyeCenterY], leftRadius,'EdgeColor','r');
%     viscircles([rightEyeCenterX, rightEyeCenterY], rightRadius,'EdgeColor', 'g');
%     pause;
    
    referenceDirection = [0 1];
    eyeDirection = [leftEyeCenterX-rightEyeCenterX, leftEyeCenterY-rightEyeCenterY];
    eyeDirection = eyeDirection / norm(eyeDirection);
    
    angle = acos(dot(referenceDirection, eyeDirection)) * 180 / pi;
    angle = 90 - angle;

%     figure; imshow(face); title('face before rotated'); pause;
    
    marker=zeros(size(face(:,:,1)));
    originalMarker = marker;
    
    face = imrotate(face, -angle);
    rgbImage = imrotate(rgbImage, -angle);
    
    marker(round(leftEyeCenterY), round(leftEyeCenterX))=1;
%     figure; imshow(marker); title('left eye marker before'); pause;
    marker_rot = imrotate(marker, -angle, 'bicubic');
    [leftEyeCenterY, leftEyeCenterX]=find(marker_rot > 0);
    leftEyeCenterY = leftEyeCenterY(1);
    leftEyeCenterX = leftEyeCenterX(1);
%     figure; imshow(marker_rot); title('left eye marker_rot'); pause;
    
    marker = originalMarker;
    marker(round(rightEyeCenterY), round(rightEyeCenterX))=1;
%     figure; imshow(marker); title('right eye marker before'); pause;
    marker_rot = imrotate(marker, -angle, 'bicubic');
%     figure; imshow(marker_rot); title('right eye marker after'); pause;
    [rightEyeCenterY, rightEyeCenterX]=find(marker_rot > 0);
    rightEyeCenterY = rightEyeCenterY(1);
    rightEyeCenterX = rightEyeCenterX(1);
%     figure; imshow(marker_rot); title('right eye marker_rot'); pause;
    
    faceMask = imrotate(faceMask, -angle);
    
%     figure; imshow(face); title('face after rotated'); pause;
%     viscircles([leftEyeCenterX, leftEyeCenterY], leftRadius,'EdgeColor','r');
%     viscircles([rightEyeCenterX, rightEyeCenterY], rightRadius,'EdgeColor', 'g');
%     pause;

end