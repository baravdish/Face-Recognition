function mask = getAttenuationMask (maskSize, radius, centerValue, edgeValue)
%
% function mask = getAttenuationMask (maskSize, radius, centerValue, edgeValue)
%

if nargin==2
  centerValue = 1;
  edgeValue = 0;
end

width = maskSize(1);
height = maskSize(2);

mx = width / 2;
my = height / 2;

mask=zeros(maskSize);

for i=1:width
  for j=1:height
      d = sqrt( (i-mx)^2 + (j-my)^2 );
      if (d >= radius)
        d = edgeValue;
      else
        d = (centerValue * (1 - (d / radius))) + (edgeValue * (d / radius));
      end

      mask(i, j) = d;
  end
end