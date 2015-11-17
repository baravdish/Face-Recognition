function [output] =  color_normalization(img)
% inp_image = imread(FILENAME); % read input image 
% http://www.inf.u-szeged.hu/ssip/2005/lectures/sergyan.pdf
inp_image = img;
[m,n,d]=size(inp_image);     % get size of input image 
f=double(inp_image);      % double needed for computations 
M=zeros(m*n,3); 
z=1; 
mv=mean(mean(f));    % a vector containing the mean r,g and b value 
v1=[mv(1),mv(2),mv(3)];    % means in red, green and blue


for i=1:m 
    for j=1:n 
        v=[f(i,j,1),f(i,j,2),f(i,j,3)];  % image pixel at i,j 
        M(z,:) = v - v1;    % image normed to mean zero 
        z = z + 1; 
    end
end
C = cov(M); % covariance computed using Matlab cov function


%find eigenvalues and eigenvectors of C. 
[V,D]=eig(C); % computes the eigenvectors(V) and eigenvalues (diagonal elements of D) of the color cluster C 
%get the max. eigenvalue meig and the corresponding eigenvector ev0.
meig = max(max(D)); % computes the maximum eigenvalue of C. Could also be norm(C) 
if(meig==D(1,1)) ev0=V(:,1); end 
if(meig==D(2,2)) ev0=V(:,2); end 
if(meig==D(3,3)) ev0=V(:,3); end % selects the eigenvector belonging to the greatest eigenvalue


Idmat =eye(3);    % identity matrix of dimension 3 
wbaxis=[1;1;1]/sqrt(3); % unit vector pointing from origin along the main diagonal 
nvec = cross(ev0,wbaxis); % rotation axis , cross(A,B)=A×B 
cosphi = dot(ev0,wbaxis) % dot product, i.e. sum((ev0.*wbaxis)) 
sinphi = norm(nvec); %  sinphi is the length of the cross product of two unit vectors 
%normalized rotation axis. 
nvec = nvec/sinphi; % normalize nvec


if(cosphi>0.99) 
    f=uint8(f); 
    output = f;
    return;
%     imwrite(f,OUTPUT); %in this case we dont normalize, output is input etc. 
else % we normalize 
    n3 = nvec(3); 
    n2 = nvec(2);
    n1 = nvec(1);   
% remember: this is a unit vector along the rotation axis
U = [[ 0  -n3  n2]; [ n3   0  -n1]; [ -n2 n1  0]]; 
U2 = U*U; 
Rphi = Idmat + (U*sinphi) + (U2*(1-cosphi));


n0   = [0 0 0]'; 
n255 = [255 255 255]'; 
for i=1:m 
    for j=1:n 
        s(1)= f(i,j,1)-mv(1); % compute vector s of normalized image at i,j 
        s(2)= f(i,j,2)-mv(2); 
        s(3)= f(i,j,3)-mv(3); 
        t = Rphi*s' ; % s transposed, as s is row vector, then rotated 
        tt = floor(t + [128 128 128]'); % shift to middle of cube and make it integer
        tt = max(tt,n0);   % handling underflow 
        tt = min(tt,n255); % handling overflow
        g(i,j,:)=tt; 
    end 
end
g=uint8(g); 
output = g;
end % end of normalization