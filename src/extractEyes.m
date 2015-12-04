function [leftEyeCenterX, leftEyeCenterY, leftEyeRadius, ...
          rightEyeCenterX, rightEyeCenterY, rightEyeRadius] = extractEyes(finalEyeMap)
      
    [centersBright, radiiBright] = imfindcircles(finalEyeMap, [4 20], ...
                                                 'ObjectPolarity', 'bright', ...
                                                 'Method', 'TwoStage', ...
                                                 'sensitivity', 0.99);
             
%     figure; imshow(I); title('I'); pause;
%     viscircles(centersDark(1:2,:), radiiDark(1:2),'EdgeColor','r');
%     viscircles(centersBright(1:2,:), radiiBright(1:2),'EdgeColor', 'b');
    
%     viscircles(centersDark(1:10,:), radiiDark(1:10), 'EdgeColor', 'r');
%     viscircles(centersDark(1:2,:), radiiDark(1:2),'EdgeColor', 'g');
%     viscircles(centersBright(1:10,:), radiiBright(1:10),'EdgeColor', 'b');

% pause;

    eyeMapSize = size(finalEyeMap);
    imageSizeX = eyeMapSize(1,1,1);
    imageSizeY = eyeMapSize(1,2,1);
    [cols rows] = meshgrid(1:imageSizeY, 1:imageSizeX);

    centerX1 = centersBright(1,1,1);
    centerY1 = centersBright(1,2,1);
    radius1 = radiiBright(1,1);
%     irisMask = ((rows - centerY1).^2 + (cols - centerX1).^2) <= radius1.^2;
    
    centerX2 = centersBright(2,1,1);
    centerY2 = centersBright(2,2,1);
    radius2 = radiiBright(2,1);
    
    for n=2 : size(centersBright, 1)
        centerX2 = centersBright(n,1,1);
        centerY2 = centersBright(n,2,1) ; 
        radius2 = radiiBright(n,1);
        
        distance = sqrt((centerX2-centerX1)^2 + (centerY2-centerY1)^2);
        
        if distance > (radius1 + radius2)*4
            n = size(centersBright, 1);
            break
        end
    end
         
%     irisMask = irisMask | ((rows - centerY2).^2 + (cols - centerX2).^2) <= radius2.^2;
%     irisMaskRep = repmat(irisMask, [1,1,3]);
%     irisMap = face.*uint8(irisMaskRep);
%     figure; imshow(irisMap); title('irisMap'); pause;
    
    leftEyeCenterX = centerX1;
    leftEyeCenterY = centerY1;
    leftEyeRadius = radius1;
    
    rightEyeCenterX = centerX2;
    rightEyeCenterY = centerY2;
    rightEyeRadius = radius2;
    if centerX1 > centerX2 
        leftEyeCenterX = centerX2;
        leftEyeCenterY = centerY2;
        leftEyeRadius = radius2;
        
        rightEyeCenterX = centerX1;
        rightEyeCenterY = centerY1;  
        rightEyeRadius = radius1;
    end
    
%     figure; imshow(face); title('face');
%     viscircles([leftEyeCenterX, leftEyeCenterY], leftRadius,'EdgeColor','r');
%     viscircles([rightEyeCenterX, rightEyeCenterY], rightRadius,'EdgeColor', 'g');
% %     viscircles([mouthX, mouthY], rightRadius,'EdgeColor', 'b');
%     pause;

end