function robustAWB(originalImg, org2)

	figure; imshow(originalImg); title('originalImg');
    
    % Convert to YUV color space
    yuvImg = rgb2yuv(originalImg);
    
%     figure; imshow(yuvImg); title('yuvImg');
    
    sRgb2XyxTransform = [0.4124564  0.3575761  0.1804375; ...
                         0.2126729  0.7151522  0.0721750; ...
                         0.0193339  0.1191920  0.9503041];
           
    rgb2yuvTransform = [0.299 0.587 0.114; ...
                        -0.299 -0.587 0.886; ...
                        0.701 -0.587 -0.114];
    
    invRgb2yuvTransform = inv(rgb2yuvTransform);
    
    % Set T higher for a larger threshold of deviation to consider off-gray
    T = 0.1; 
    F = ( abs(yuvImg(2,:)) + abs(yuvImg(3,:)) )./yuvImg(1,:);
    grays = yuvImg(:,F<T);
    U_bar = mean(grays(2,:));           %estimate of the U of the illuminant
    V_bar = mean(grays(3,:));           %estimate of the V of the illuminant
%     size([100 ,11 ,11])
%     size(invRgb2yuvTransform)

    rgbEst = invRgb2yuvTransform*[100;U_bar;V_bar]; %convert the average gray from YUV to RGB

    xyEst = xyz2xy(sRgb2XyxTransform*rgbEst);   %calculate xy chromaticity
    
    xyzEst = xyz2xy([xyEst,100]);         %normalize Y to 100 so D65 luminance comparable

    
    % gain adjustment parameters, can probably be optimized via more control system analysis
    u = .01; %gain step size
    a = .8; %double step threshold
    b = .001; %convergence threshold
    

    theSize = size(originalImg)
    imRGB_orig = reshape(originalImg, theSize(1) * theSize(2), 3);
    size(imRGB_orig)
    imRGB = imRGB_orig;
    gain = [1 1 1];
    
    maxIter = 100;
    catType = 'cat02';
    xyz_D65 = [95.04; 100; 108.88]; %http://en.wikipedia.org/wiki/D65, normalized Y = 100
    
    U_avg = [];
    V_avg = [];
    totGray = [];
    for iter = 1:maxIter
        
        im = double(imRGB)*rgb2yuvTransform; %convert to YUV
        % find gray chromaticity
        % (|U|+|V|)/Y
        F = ( abs(im(2,:)) + abs(im(3,:)) )./im(1,:);

    %     figure
    %     imshow(cbunshape(imRGB/255.*repmat(F<T,3,1),size(im_orig)))
    %     % imshow(cbunshape(inv_xfm)*im+[F>T;zeros(2,length(F))],size(im_orig))) %rubylith overlays
    %     title('Gray Pixels Found by Robust AWB'),xlabel('Pixels in their original color are considered gray')

        totGray = [totGray sum(F<T)];
        if totGray(end) == 0
            disp('No valid gray pixels found.')
            break
        end
        grays = im(:,F<T);
        U_bar = mean(grays(2,:));
        V_bar = mean(grays(3,:));

        U_avg = [U_avg U_bar];
        V_avg = [V_avg V_bar];


        if max(abs([U_bar V_bar])) < b
            disp(['Converged. U_bar and V_bar < ' num2str(b) ' in magnitude.'])
            break
        elseif iter >= 2 && norm([U_avg(end)-U_avg(end-1) V_avg(end)-V_avg(end-1)]) < 10^-6
            disp(['U_bar and V_bar are no longer improving.'])
            break
        end
        rgbEst = invRgb2yuvTransform*[100;U_bar;V_bar]; %convert the average gray from YUV to RGB
        xyEst = xyz2xy(sRgb2XyxTransform*rgbEst); %calculate xy chromaticity
        xyzEst = xy2xyz(xyEst,100); %normalize Y to 100 so D65 luminance comparable
        imRGB = double(imRGB)*cbCAT(xyzEst, xyz_D65, catType);


%         if abs(U_bar) > abs(V_bar) % U > V; blue needs adjustment
%             err = U_bar;
%             ch = 3; %blue channel
%         else
%             err = V_bar;
%             ch = 1; %red channel
%         end
%         if abs(err) >= a
%             delta = 2*sign(err)*u; %accelerate gain adjustment if far off
%         elseif abs(err) < b %converged
%             delta = 0;
%             disp(['Converged. U_bar and V_bar < ' num2str(b) ' in magnitude.'])
%             break
%         else
%             delta = err*u;
%         end
% 
%         gain(ch) = gain(ch)-delta; %negative fdbk loop
%         imRGB = double(imRGB_orig) * diag(gain);
    end
    imRGB = reshape(imRGB, theSize(1), theSize(2), 3);
    
    figure; imshow(imRGB); title('imRGB');
    figure; imshowpair(originalImg,imRGB,'blend','Scaling','joint'); title('blend1');
%     figure; imshow(imabsdiff(originalImg,imRGB),[]); title('blend2');
    
%     M = rand(numRows, numColumns, 3);
%     T = rand(3, 3);
%     P = reshape(M, numRows * numColumns, 3) * T';
%     Result = reshape(P, numRows, numColumns, 3);
    
    figure
    plot(totGray)
    title('rAWB: Total Gray Pixels vs Iterations')
    print(gcf,'-dpng',[outFile '-fig1'])
%     
%     figure
%     hold on
%     plot(U_avg,'b')
%     plot(V_avg,'r')
%     hold off
%     title('rAWB: Mean Chromaticity vs Iterations')
%     legend('U','V')
%     print(gcf,'-dpng',[outFile '-fig2'])
        
%     figure
%     imshow(cbunshape(imRGB/255.*repmat(F<T,3,1),size(im_orig)))
%     % imshow(cbunshape(inv_xfm)*im+[F>T;zeros(2,length(F))],size(im_orig))) %rubylith overlays
%     title('Gray Pixels Found by Robust AWB'),xlabel('Pixels in their original color are considered gray')
%     print(gcf,'-dpng',[outFile '-fig3'])

end

% 
% 
% function robustAWB(filename,outFile,option,catType,T,maxIter,plots)
% %robustAWB(filename,outFile,option,catType,T,maxIter,plot)
% % Performs robust auto white-balancing by estimating gray pixels based on
% % their deviation in YUV space then applying an iterative correction via
% % CAT or directly adjusting the R and B channels.
% % Set option = 'RB gain' or 'cat' for the correction method.
% % Set T higher for a larger threshold of deviation to consider off-gray.
% % Set plot = 0 or 1 to turn diagnostic plots on or off.
% 
% tutorialinit
% 
% %% import image
% im_orig = cbimread(filename);
% % figure
% % imshow(im_orig)
% % title('Original Image')
% 
% %% various constants
% xyz_D65 = [95.04; 100; 108.88]; %http://en.wikipedia.org/wiki/D65, normalized Y = 100
% 
% sRGBtoXYZ =   [0.4124564  0.3575761  0.1804375; ...
%                0.2126729  0.7151522  0.0721750; ...
%                0.0193339  0.1191920  0.9503041];
% 
% %% robust auto white balance - YUV grays
% % option = 'RB gain'; % adjustments are done as a gain to only R or B channelas descripted in Huo
% % option = 'cat';
% % catType = 'CAT02';
% 
% % maxIter = 1000;
% 
% % T = 0.3; %this should be tweakable slider
% % T = 0.1;
% % gain adjustment parameters, can probably be optimized via more control system analysis
% u = .01; %gain step size
% a = .8; %double step threshold
% b = .001; %convergence threshold
% 
% % rgb to yuv
% xfm =   [0.299 0.587 0.114; ...
%         -0.299 -0.587 0.886; ...
%         0.701 -0.587 -0.114];
%     
% inv_xfm = inv(xfm);
%     
% imRGB_orig = cbreshape(im_orig)*255;
% imRGB = imRGB_orig;
% gain = [1 1 1];
% 
% U_avg = [];
% V_avg = [];
% totGray = [];
% for iter = 1:maxIter
%     im = xfm*imRGB; %convert to YUV
%     % find gray chromaticity
%     % (|U|+|V|)/Y
%     F = ( abs(im(2,:)) + abs(im(3,:)) )./im(1,:);
%     
% %     figure
% %     imshow(cbunshape(imRGB/255.*repmat(F<T,3,1),size(im_orig)))
% %     % imshow(cbunshape(inv_xfm)*im+[F>T;zeros(2,length(F))],size(im_orig))) %rubylith overlays
% %     title('Gray Pixels Found by Robust AWB'),xlabel('Pixels in their original color are considered gray')
%     
%     totGray = [totGray sum(F<T)];
%     if totGray(end) == 0
%         disp('No valid gray pixels found.')
%         break
%     end
%     grays = im(:,F<T);
%     U_bar = mean(grays(2,:));
%     V_bar = mean(grays(3,:));
%     
%     U_avg = [U_avg U_bar];
%     V_avg = [V_avg V_bar];
%     
%     if strcmpi(option,'cat')
%         if max(abs([U_bar V_bar])) < b
%             disp(['Converged. U_bar and V_bar < ' num2str(b) ' in magnitude.'])
%             break
%         elseif iter >= 2 && norm([U_avg(end)-U_avg(end-1) V_avg(end)-V_avg(end-1)]) < 10^-6
%             disp(['U_bar and V_bar are no longer improving.'])
%             break
%         end
%         rgbEst = inv_xfm*[100;U_bar;V_bar]; %convert the average gray from YUV to RGB
%         xyEst = XYZ2xy(sRGBtoXYZ*rgbEst); %calculate xy chromaticity
%         xyzEst = xy2XYZ(xyEst,100); %normalize Y to 100 so D65 luminance comparable
%         imRGB = cbCAT(xyzEst,xyz_D65,catType)*imRGB;
%     else
%         if abs(U_bar) > abs(V_bar) % U > V; blue needs adjustment
%             err = U_bar;
%             ch = 3; %blue channel
%         else
%             err = V_bar;
%             ch = 1; %red channel
%         end
%         if abs(err) >= a
%             delta = 2*sign(err)*u; %accelerate gain adjustment if far off
%         elseif abs(err) < b %converged
%             delta = 0;
%             disp(['Converged. U_bar and V_bar < ' num2str(b) ' in magnitude.'])
%             break
%         else
%             delta = err*u;
%         end
% %         [err,delta]
%         gain(ch) = gain(ch)-delta; %negative fdbk loop
%         imRGB = diag(gain)*imRGB_orig;
%     end
% end
% 
% imwrite(cbunshape(imRGB,size(im_orig))/255,outFile,'png');
% % figure
% % imshow(cbunshape(imRGB,size(im_orig))/255)
% % title('Robust AWB Corrected')
% if plots
%     length(U_avg)
%     figure
%     plot(totGray)
%     title('rAWB: Total Gray Pixels vs Iterations')
%     print(gcf,'-dpng',[outFile '-fig1'])
%     
%     figure
%     hold on
%     plot(U_avg,'b')
%     plot(V_avg,'r')
%     hold off
%     title('rAWB: Mean Chromaticity vs Iterations')
%     legend('U','V')
%     print(gcf,'-dpng',[outFile '-fig2'])
%         
%     figure
%     imshow(cbunshape(imRGB/255.*repmat(F<T,3,1),size(im_orig)))
%     % imshow(cbunshape(inv_xfm)*im+[F>T;zeros(2,length(F))],size(im_orig))) %rubylith overlays
%     title('Gray Pixels Found by Robust AWB'),xlabel('Pixels in their original color are considered gray')
%     print(gcf,'-dpng',[outFile '-fig3'])
% end