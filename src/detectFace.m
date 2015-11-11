%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: Image in which to detect face.
%
% output: Image containing only the face.
%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = detectFace(input)
  
    imgYCbCr = rgb2ycbcr(input);
%     imgHSV = rgb2hsv(input);
    
    Y = imgYCbCr(:,:,1);
    CB = imgYCbCr(:,:,2);
    CR = imgYCbCr(:,:,3);
    
%     figure
%     fontSize = 12;

%   subplot(4, 3, 1);
% 	imshow(Y);
% 	title('Y', 'FontSize', fontSize);
% 	subplot(4, 3, 2);
% 	imshow(CB);
% 	title('CB', 'FontSize', fontSize);
% 	subplot(4, 3, 3);
% 	imshow(CR);
% 	title('CR', 'FontSize', fontSize);
%     
%     subplot(4, 3, 4);
% 	plot(imhist(Y));
% 	title('hist Y', 'FontSize', fontSize);
% 	subplot(4, 3, 5);
% 	plot(imhist(CB));
% 	title('hist CB', 'FontSize', fontSize);
% 	subplot(4, 3, 6);
% 	plot(imhist(CR));
% 	title('hist CR', 'FontSize', fontSize);
%     
%     subplot(4, 3, 7);
% 	imshow(H);
% 	title('H', 'FontSize', fontSize);
% 	subplot(4, 3, 8);
% 	imshow(S);
% 	title('S', 'FontSize', fontSize);
% 	subplot(4, 3, 9);
% 	imshow(V);
% 	title('V', 'FontSize', fontSize);
%     
%     subplot(4, 3, 10);
% 	plot(imhist(H));
% 	title('hist H', 'FontSize', fontSize);
% 	subplot(4, 3, 11);
% 	plot(imhist(S));
% 	title('hist S', 'FontSize', fontSize);
% 	subplot(4, 3, 12);
% 	plot(imhist(V));
% 	title('hist V', 'FontSize', fontSize);
    
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

%     figure
%     title('faceMask' + int2str(i))
%     imshow(faceMask);
    
    faceMask = repmat(faceMask, [1,1,3]);
    output = img.*uint8(faceMask);

end
