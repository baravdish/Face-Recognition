%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: Image in which to detect face.
%
% output: Image containing only the face.
%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = detectFace(input)
  
    imgYCbCr = rgb2ycbcr(input);
    
    Y = imgYCbCr(:,:,1);
    CB = imgYCbCr(:,:,2);
    CR = imgYCbCr(:,:,3);

    faceMask = and(and(CB > 105, CB < 135), ... 
    and(CR > 140, CR < 165));

    for n = 1:6
        erodeKernel = strel('disk', n);
        dilateKernel = strel('disk', 1+n);
        faceMask = imerode(faceMask, erodeKernel);
        faceMask = imdilate(faceMask, dilateKernel);
    end
   
    erodeKernel = strel('disk', 25);
    dilateKernel = strel('disk', 20);    
    faceMask = imdilate(faceMask, dilateKernel);
    faceMask = imerode(faceMask, erodeKernel);

    erodeKernel = strel('disk', 10);
    dilateKernel = strel('disk', 10); 
    faceMask = imerode(faceMask, erodeKernel);
    faceMask = imdilate(faceMask, dilateKernel);
    
    faceMask = repmat(faceMask, [1,1,3]);
    face = input.*uint8(faceMask);
    
    [row, col] = find(face(:,:,1) ~= 0);
    
    minCol = min(col);
    maxCol = max(col);

    minRow = min(row);
    maxRow = max(row);
    cropImg = input(minRow:maxRow, minCol:maxCol, :);

    % EyeMap & MouthMap
    faceYCbCr = rgb2ycbcr(cropImg);
     
    faceY = double(faceYCbCr(:,:,1));
    faceCB = double(faceYCbCr(:,:,2));
    faceCR = double(faceYCbCr(:,:,3)); 

    sqFaceCB = faceCB.^2;
    sqFaceCR = faceCR.^2;
    cHat = 255 - faceCR;
    sqCHat = cHat.^2;
    divFace = faceCB./faceCR;
    
    % Now it's normalized to [0,1]
    % If we need [0,255], just multiply with 255.
    normFaceCB = (sqFaceCB)/(max(sqFaceCB(:)));
    normFaceCR = (sqFaceCR)/(max(sqFaceCR(:)));
    
    normSqCHat = (sqCHat)/(max(sqCHat(:)));
    normDivFace = (divFace)/(max(divFace(:)));
    
    eyeMapC = 1/3*(normFaceCB + normSqCHat + normDivFace);
    eyeMapCEnhanc = imadjust(eyeMapC);
    
    diff = abs(normFaceCB - normDivFace);
    mouthMap = and(diff, normFaceCR);
    
    figure
    imshow(eyeMapCEnhanc);
    output = cropImg;

end
