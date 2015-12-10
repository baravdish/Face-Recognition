function [estimatedSkinMask] = extractSkinColorFromCbCr(Cb, Cr)
   
    % Values inspired by:
    % https://web.stanford.edu/class/ee368/Project_03/Project/reports/ee368group15.pdf
    
    estimatedSkinMask = and(and(Cb > 95, Cb < 145), ... 
                            and(Cr > 132, Cr < 165));
    
end