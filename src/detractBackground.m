function foregroundMask = detractBackground(img)

X = imadjust(rgb2gray(img));
%     for n = 1:10
      figure;
%       level = graythresh(IDX);
%       im2bw(IDX, level)
%       masked = X > ( 0.1*255 * (n-1) );
%       masked = X < masked + 0.1*255;
%     X = imfilter(X, fspecial('gaussian',2));

        overSaturated = rgb2gray(img);
%     imshow(overSaturated); title('Over saturation'); pause;
%     imshow(imadjust(overSaturated)); title('imadjust(overSaturated)'); pause;
    imshow(imfilter(imadjust(overSaturated), fspecial('laplacian'))); title('filter imadjust(overSaturated)'); pause;
    overSaturated = imfilter(overSaturated, fspecial('laplacian'));
%     imshow(overSaturated); title('filter Over saturation'); pause;
    overSaturated = imdilate(overSaturated, strel('disk', 2));
    imshow(overSaturated); title('dilate Over saturation'); pause;
    overSaturated = ~overSaturated;
    imshow(overSaturated); title('filter Over saturation'); pause;
    imshow(X); title('X before saturation'); pause;
%     overSaturatedRep = repmat(~overSaturated, [1,1,3]);
    X = X.*uint8(~overSaturated);
    imshow(X); title('X after saturation'); pause;
    
    masked = edge(X, 'canny');
    
%     mask = zeros(size(X));
%     mask(25:end-25,25:end-25) = 1;
%     bw = activecontour(rgb2gray(img),mask, 300);

%     masked = imgradient(X) / 255;
%     imshow(masked);
%     title('imgradient');
%     pause;
%     
%     hsvImg = rgb2hsv(img);
%     imshow(hsvImg);
%     title('hsvImg');
%     pause;
    
%  
    
    
      imshow(masked);
%       title(['n = ' int2str(n)]);
      pause;
      
      masked = imdilate(masked, strel('disk', 2));
      blobs = ~masked;
      blobs = imerode(blobs, strel('disk', 2));
      
      imshow(blobs);
      pause;
      theSize = size(img,1)*size(img,2)*0.05
%       overSaturated = 1000
      blobs = bwareaopen(blobs, floor(theSize));
      imshow(blobs);
      pause;
      

      
      lapImg = X;
%       lapImg = imfilter(lapImg, fspecial('gaussian',10));
        lapImg = medfilt2(lapImg);

      lapImg = imfilter(lapImg, fspecial('log'));
%       lapImg = lapImg.*lapImg;
%       lapImg = imadjust(lapImg);
      background = zeros(size(lapImg));
      imshow(lapImg); title('lapImg'); pause;
      
      [labeledImage, numberOfBlobs] = bwlabel(blobs);
       blobMeasurements = regionprops(labeledImage, 'area');
       % Get all the areas
       allAreas = [blobMeasurements.Area];
	  [sortedAreas, sortIndexes] = sort(allAreas, 'descend');

      
      for n2 = 1:length(sortIndexes)
          largest = ismember(labeledImage, sortIndexes(n2));
%           blobs = blobs - largest;
          imshow(largest); pause;

          
          foreground = largest;
          size(foreground)
          length(nonzeros(foreground(:)))
          
          rega = lapImg.*uint8(foreground);
          percentage = (100 * sum(rega(:))) / (length(nonzeros(foreground(:)))*255)
          imshow(rega); title('rega'); pause;
          if percentage < 0.8
              background = background | foreground;
          end
      end
      imshow(background); title('background'); pause;
        backgroundRep = repmat(~background, [1,1,3]);
    foreground = img.*uint8(backgroundRep);
      imshow(foreground); title('foreground'); pause;
      
      foregroundMask = ~background;
end