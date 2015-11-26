%%%%%%%%%%%%%%%%%%%%%%%%%
% input: RGB-Image to correct.

% output: Corrected image.
%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = colorCorrection(img)
    
    R = img(:,:,1);
    G = img(:,:,2);
    B = img(:,:,3);

    rVec = R(:);
    gVec = G(:);
    bVec = B(:);

    sortedImg(:,1) = sort(rVec, 'descend');
    sortedImg(:,2) = sort(gVec, 'descend'); 
    sortedImg(:,3) = sort(bVec, 'descend');

    nPixels = length(rVec);
    whitePercentage = 0.001;
    nFindWhite = round(nPixels*whitePercentage);
    blackPercentage = 0.015;
    nFindBlack = round(nPixels*blackPercentage);

    brighMean = [mean(sortedImg((1:nFindWhite),1)), ...
                 mean(sortedImg((1:nFindWhite),2)), ...
                 mean(sortedImg((1:nFindWhite),3))];

    darkMean = [mean(sortedImg((end-nFindBlack:end),1)), ...
                mean(sortedImg((end-nFindBlack:end),2)), ...
                mean(sortedImg((end-nFindBlack:end),3))];

    [rows, columns] = size(R);
    R = R - (ceil(darkMean(1,1))*uint8(ones(rows, columns)));
    G = G - (ceil(darkMean(1,2))*uint8(ones(rows, columns)));
    B = B - (ceil(darkMean(1,3))*uint8(ones(rows, columns)));
     
    % figure; imshow(img); title('before'); pause;   
    % output = cat(3, R, G, B);
    % imshow(output); title('after'); pause;   
    % return;

    scaleFactors = 255./brighMean;
    
    rNorm = R*scaleFactors(1);
    binR = rNorm > 255;
    rNorm(binR) = 255;

    gNorm = G*scaleFactors(2);
    binG = gNorm > 255;
    gNorm(binG) = 255;

    bNorm = B*scaleFactors(3);
    binB = bNorm > 255;
    bNorm(binB) = 255;
    output = cat(3, rNorm, gNorm, bNorm);

    % figure; imshow(originalImg); title('original'); pause;
    % figure; imshow(output); title('WhiteBalanced'); pause;

end
