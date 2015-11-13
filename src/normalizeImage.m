%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: Image to be normalized.
%
% output: Normalized image.
%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = normalizeImage(input)
    faceYCbCr = rgb2ycbcr(input);
    
    faceY = double(faceYCbCr(:,:,1));
    faceCB = double(faceYCbCr(:,:,2));
    faceCR = double(faceYCbCr(:,:,3)); 

    sqFaceCB = faceCB.^2;
    sqFaceCR = faceCR.^2;
    cHat = 255 - faceCR;
    sqCHat = cHat.^2;
    divFace = faceCB./faceCR;
    
    % Normalized to [0,1]
    normFaceCB = (sqFaceCB)/(max(sqFaceCB(:)));
    normFaceCR = (sqFaceCR)/(max(sqFaceCR(:)));
    
    normSqCHat = (sqCHat)/(max(sqCHat(:)));
    normDivFace = (divFace)/(max(divFace(:)));
    
    eyeMapC = 1/3*(normFaceCB + normSqCHat + normDivFace);
    eyeMapCEnhanc = imadjust(eyeMapC);
    figure
    imshow(eyeMapCEnhanc);
    
    
    
%     diff = abs(normFaceCB - normDivFace);
%     mouthMap = and(diff, normFaceCR);
    output = eyeMapCEnhanc;
    
end
