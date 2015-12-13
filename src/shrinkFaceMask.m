function faceMask = shrinkFaceMask(faceMask) 
    faceMaskCopy = faceMask;
    minCol = [];
    maxCol = [];
    minRow = [];
    maxRow = [];
    eraseRadius = 80;
    maxIterations = 100;
    while( (isempty(minCol) || isempty(maxCol) || ... 
           isempty(minRow) || isempty(maxRow)) && maxIterations > 0)
        faceMaskCopy = faceMask;
        minCol = [];
        maxCol = [];
        minRow = [];
        maxRow = [];
        faceMaskCopy = imerode(faceMaskCopy, strel('disk', eraseRadius));
        faceMaskCopy = imdilate(faceMaskCopy, strel('disk', eraseRadius));
        [row, col] = find(faceMaskCopy(:,:,1) ~= 0);
        minCol = min(col);
        maxCol = max(col);
        minRow = min(row);
        maxRow = max(row);
        eraseRadius = max(eraseRadius - 10, 1);
        maxIterations = maxIterations - 1;
    end
    faceMask = faceMaskCopy;
end