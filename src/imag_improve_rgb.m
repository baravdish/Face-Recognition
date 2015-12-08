%%%%%%%%%%%%%%%RGB normalisation%%%%%%%%%%%%%%%%%%%%%%
%its cascaded implementain of 1 section of paper "A FAST SKIN REGION DETECTOR" by
%Phil Chen, Dr.Christos Greecos
%and
%section 2.1 of paper "Simple and accurate face detection in color images" by
%YUI TING PAI et al
% Coding by Madhava.S.Bhat, Dept of Electronics an communication
%Dr.Ambedkar Institute of Technology, Bangalore
%madhava.s@dr-ait.org
function Cfinal = imag_improve_rgb(IMG)
figure,imshow(IMG)
title('original')
R=double(IMG(:,:,1));
G=double(IMG(:,:,2));
B=double(IMG(:,:,3));
[H,W]=size(R);

[srow,scol]=find(R==0 & G==0 & B==0);
% if(isempty(srow) && isempty(scol))
    minR=min(min(R));
    minG=min(min(G));
    minB=min(min(B));
% end
R=R-minR;
G=G-minG;
B=B-minB;

S=zeros(H,W);
[srow,scol]=find(R==0 & G==0 & B==0);
[sm,sn]=size(srow);

for i=1:sm
     S(srow(i),scol(i))=1;
end
mstd=sum(sum(S));
Nstd=(H*W)-mstd;

Cst=0;
Cst=double(Cst);
for i=1:H
    for j=1:W
         a=R(i,j);
         b=R(i,j);
         
            if(B(i,j)<a)
               a=B(i,j);
            else
               b=B(i,j);
            end
         
             if(G(i,j)<a)
                a=G(i,j);
             else
                b=G(i,j);
             end
         
         Cst=a+b+Cst;
         
    end
end
%%%%sum of black pixels%%%%%%%%%%%
blacksumR=0;
blacksumG=0;
blacksumB=0;
for i=1:sm
    blacksumR=blacksumR+R(srow(i),scol(i));
    blacksumG=blacksumG+G(srow(i),scol(i));
    blacksumB=blacksumR+B(srow(i),scol(i));
end
Cstd = Cst/(2*Nstd);
CavgR=sum(sum(R))./(H*W);
CavgB=sum(sum(B))./(H*W);
CavgG=sum(sum(G))./(H*W);
Rsc=Cstd./CavgR;
Gsc=Cstd./CavgG;
Bsc=Cstd/CavgB;

R=R.*Rsc;
G=G.*Gsc;
B=B.*Bsc;

C(:,:,1)=R;
C(:,:,2)=G;
C(:,:,3)=B;
C=C/255;
YCbCr=rgb2ycbcr(C);
Y=YCbCr(:,:,1);

figure,imshow(C)
title('aft 1st stage of compensation')

%normalize Y
minY=min(min(Y));
maxY=max(max(Y));
Y=255.0*(Y-minY)./(maxY-minY);
YEye=Y;
Yavg=sum(sum(Y))/(W*H);

T=1;
if (Yavg<64)
    T=1.4;
elseif (Yavg>192)
    T=0.6;
end

if (T~=1)
    RI=R.^T;
    GI=G.^T;
else
    RI=R;
    GI=G;
end


Cfinal(:,:,1)=uint8(RI);
Cfinal(:,:,2)=uint8(GI);
Cfinal(:,:,3)=uint8(B);
figure,imshow(Cfinal)
title('Light intensity compensated')

% YCbCr=rgb2ycbcr(Cnew);