 function maxellipse(Image)

    ImageArea=nnz(Image);

    [m,n]=size(Image);
    cxy=[m,n]/2+.5;

    [x,y]=ndgrid((1:m)-cxy(1),(1:n)-cxy(2));
     xygrid=[x(:),y(:)].';

     ImgAdjusted=single(Image*100-99);
     fun=@(p) objective(p,ImgAdjusted,xygrid, ImageArea); %tuning params

               %%Initial Guess
              S=regionprops(Image,'Orientation',...
                           'MajorAxisLength','MinorAxisLength' );

              p0(1)=S.MajorAxisLength/3;
              p0(2)=S.MinorAxisLength/3;
              p0(3)=S.Orientation;     

     p=fminsearch(fun,p0(:));

     [~,report]=fun(p);

     imagesc(report.ellipsemask+Image);

function  [cost,report]=objective(p,Image,xygrid, ImageArea)

   a=abs(p(1)); %ellipse length a 
   b=abs(p(2)); %ellipse length b
   theta=p(3);

   R=@(t)[cosd(t), -sind(t); sind(t) cosd(t) ];
    R=R(-theta);

   D=diag(1./([a,b]/2).^2);

   strel =  sum( ((R*D*R.')*xygrid).*xygrid, 1 ) <=1; %mask
    strel=reshape(single(strel), size(Image));

    map=conv2(Image,strel,'same');
    overlap=max(map(:));

    area=min(pi*a*b, ImageArea );

    cost=-(area+overlap);

    if nargout>1
        Map=conv2(double(Image>0),strel,'same')>=nnz(strel);
        idx=find(Map);
        Map(idx(2:end))=0;
        ImageEroded=Image & ~(conv2(double(Map),strel,'same')>0);
        BW=Image& ~ImageEroded;
        report=regionprops(BW,'Centroid','Orientation',...
                              'MajorAxisLength','MinorAxisLength' );
        report.ellipsemask=BW;
    end