function mouthMap = extractMouthMap(faceMask, face)

    ycbcrFace = rgb2ycbcr(face);
    
    Cb = double(ycbcrFace(:,:,2));
    Cr = double(ycbcrFace(:,:,3));
    
    n = length(nonzeros(faceMask));
    
    CrSquared = Cr.^2;
    numerator = (1 / n) * sum(CrSquared(:));
    
    CrDividedByCb = Cr./Cb;
    denumerator = (1 / n) * sum(CrDividedByCb(:));
    
    nn = 0.95 * (numerator / denumerator);
    
    mouthMap = CrSquared.*(CrSquared - nn*CrDividedByCb).^2;
    mouthMap = im2double(mouthMap);
    mouthMap = 255 * ((mouthMap - min(mouthMap(:))) ...
                 ./ (max(mouthMap(:)) - min(mouthMap(:))));
    mouthMap = uint8(mouthMap);
    mouthMap = mouthMap .* uint8(faceMask);
    
end