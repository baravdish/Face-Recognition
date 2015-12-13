function nonFaceMask = extractNonFaceMask(Cb, Cr, hsvImage)
    
    S = hsvImage(:,:,2);
    V = hsvImage(:,:,3);
    
    nonSkinMask1 = Cb > Cr + 10;

    nonSkinMask2 = V < S - 0.2;
    
    nonFaceMask = or(nonSkinMask1, nonSkinMask2);
    nonFaceMask = imerode(nonFaceMask, strel('disk', 5));
    nonFaceMask = bwareaopen(nonFaceMask, 1000);
    nonFaceMask = imdilate(nonFaceMask, strel('disk', 5));

end