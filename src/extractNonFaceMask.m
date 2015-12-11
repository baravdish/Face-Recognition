function nonFaceMask = extractNonFaceMask(Cb, Cr, hsvImage)
    
    S = hsvImage(:,:,2);
    V = hsvImage(:,:,3);
    
    nonSkinMask1 = Cb > Cr + 10;
%     figure; imshow(nonSkinMask1); title('nonSkinMask1'); pause;

    nonSkinMask2 = V < S - 0.2;
%     figure; imshow(nonSkinMask2); title('nonSkinMask2'); pause;
    
    nonFaceMask = or(nonSkinMask1, nonSkinMask2);
    nonFaceMask = imerode(nonFaceMask, strel('disk', 5));
    nonFaceMask = bwareaopen(nonFaceMask, 1000);
    nonFaceMask = imdilate(nonFaceMask, strel('disk', 5));
%     figure; imshow(nonFaceMask); title('nonFaceMask'); pause;

end