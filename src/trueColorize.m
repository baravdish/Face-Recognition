function binaryImage = trueColorize(rgbImage, mask, tolerance)

    
    x = rgbImage(mask);
    num_el = numel(x);
    n = 3;
    x(numel(x) + (n - mod(numel(x), n))) = 0;

    rgb_columns = reshape(x, [], n);
%     rgb_columns = reshape(rgbImage(mask), [], 3);


    size(rgb_columns)
    
    [unique_colors, m, n] = unique(rgb_columns, 'rows');
    color_counts = accumarray(n, 1);

    fprintf('There are %d unique colors in the image.\n', ...
        size(unique_colors, 1))
   
    [max_count, idx] = max(color_counts);
    
    cR = unique_colors(idx, 1)
    cG = unique_colors(idx, 2)
    cB = unique_colors(idx, 3)

    fprintf('The color [%d %d %d] occurs %d times.\n', ...
        unique_colors(idx, 1), unique_colors(idx, 2), ...
        unique_colors(idx, 3), max_count)
    
    bw = n == idx;
    
    common = rgb2ycbcr([cR, cG, cB])
    
    cR = common(1);
    cG = common(2);
    cB = common(3);
    
    rgbImage = rgb2ycbcr(im2double(rgbImage));
    
    binaryImage = (cR - tolerance(1)) < rgbImage(:,:,1) & rgbImage(:,:,1) < (cR + tolerance(1)) & ...
                  (cG - tolerance(2)) < rgbImage(:,:,2) & rgbImage(:,:,2) < (cG + tolerance(2)) & ...
                  (cB - tolerance(3)) < rgbImage(:,:,3) & rgbImage(:,:,3) < (cB + tolerance(3)) ;
    
%     binaryImage = (cR - tolerance) < rgbImage(:,:,1) & rgbImage(:,:,1) < (cR + tolerance) & ...
%                   (cG - tolerance) < rgbImage(:,:,2) & rgbImage(:,:,2) < (cG + tolerance) & ...
%                   (cB - tolerance) < rgbImage(:,:,3) & rgbImage(:,:,3) < (cB + tolerance) ;
     
%     figure; imshow(binaryImage); title('binaryImage'); pause;
    


end