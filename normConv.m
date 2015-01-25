function R=normConv(I,sigma)

I=double(I);
[dimx,dimy]=size(I);
R=zeros(dimx,dimy);

winRadius=3*sigma-1;
winSize=2*winRadius+1;

if(winSize>dimx || winSize>dimy)
    error('Kernal is too big!');
end

hg=fspecial('gaussian',winSize,sigma);

for i=1:1:dimx
    flagx=0;
    xa=max(1, i-winRadius);
    xb=min(dimx, i+winRadius);
    if(xa~=(i-winRadius))
        kxa=2-i+winRadius;
        kxb=winSize;
        flagx=1;
    elseif(xb~=(i+winRadius))
        kxa=1;
        kxb=winSize-(i+winRadius-dimx);
        flagx=1;
    else
        kxa=1;
        kxb=winSize;
    end
    
    for j=1:1:dimy
        flagy=0;
        ya=max(1, j-winRadius);
        yb=min(dimy, j+winRadius);
        if(ya~=(j-winRadius))
            kya=2-j+winRadius;
            kyb=winSize;
            flagy=1;
        elseif(yb~=(j+winRadius))
            kya=1;
            kyb=winSize-(j+winRadius-dimy);
            flagy=1;
        else
            kya=1;
            kyb=winSize;
        end
        
         tmp=hg(kxa:1:kxb, kya:1:kyb);
        if(flagx || flagy)         
            tmp=tmp./sum(sum(tmp));
        end
        conMat = tmp.*I(xa:1:xb, ya:1:yb);
        
        R(i,j)=sum(conMat(:));
        
    end
end