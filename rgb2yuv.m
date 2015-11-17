function dst = rgb2yuv(src)

% ensure this runs with rgb images as well as rgb triples
if(length(size(src)) > 2)
    
    % rgb image ([r] [g] [b])
    r = double(src(:,:,1));
    g = double(src(:,:,2));
    b = double(src(:,:,3));
    
elseif(length(src) == 3)
    
    % rgb triplet ([r, g, b])
    r = double(src(1));
    g = double(src(2));
    b = double(src(3));
    
else
    
    % unknown input format
    error('rgb2yuv: unknown input format');
    
end
    
% convert...
y = floor(0.3*r + 0.5881*g + 0.1118*b);
u = floor(-0.15*r - 0.2941*g + 0.3882*b + 128);
v = floor(0.35*r - 0.2941*g - 0.0559*b + 128);

% ensure valid range for uint8 values
y(y > 255) = 255;
y(y < 0)   = 0;
u(u > 255) = 255;
u(u < 0)   = 0;
v(v > 255) = 255;
v(v < 0)   = 0;


% generate output
if(length(size(src)) > 2)
    
    % yuv image ([y] [u] [v])
    dst(:,:,1) = uint8(y);
    dst(:,:,2) = uint8(u);
    dst(:,:,3) = uint8(v);
    
else
    
    % yuv triplet ([y, u, v])
    dst = uint8([y, u, v]);
    
end