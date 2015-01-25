function ctl = checkCircle(ctl,opt)

cc= bwconncomp(ctl);
ep= bwmorph(ctl,'endpoint');
[dimx,dimy]=size(ctl);

for i=1:1:cc.NumObjects
    a=cc.PixelIdxList{i};
    tmp=ep(a);
    t=nnz(tmp);
    if(t~=2)
        if(t==0)
            timg=zeros(dimx,dimy);
            timg(a)=1;
            idxList=sortPixels(timg, dimx, dimy);
            [newImg,flag]=Curvature(idxList,timg,opt);
            
            sig=0;
            while(~flag)
                opt.maxCurvature = 0.9*opt.maxCurvature;
                [newImg,flag]=Curvature(idxList,timg,opt);
                sig=sig+1;
                if(sig>10)
                    disp('inf loop');
                    keyboard;
                end
            end
            
            ctl(a)=0;
            ctl(newImg)=1;
        else
            disp('error in number of ep');
            keyboard
        end
    end
end