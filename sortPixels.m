function idxList=sortPixels(R, dimx, dimy)

R0=R;
nbr=[-1,-1;-1,0;-1,1;0,-1;0,1;1,-1;1,0;1,1];

%%%% check whether R is one-pixel wide %%%%

%%%% prune the boundary first %%%%
ep=bwmorph(R,'endpoint');
iter=0;
while(nnz(ep)>0)
    if(iter>10)
        disp('error in boundary');
        keyboard;
    end
    R=bwmorph(R,'spur');
    ep=bwmorph(R,'endpoint');
    iter=iter+1;
end
R=bwmorph(R,'thin',1);
%%%%%%%%%%%%%%%%%%%

tmp=[1,1,1;1,0,1;1,1,1];
numPixel = nnz(R);
labels=conv2(double(R),tmp,'same');
labels(~R)=0;
if(any(labels(:)==1) || any(labels(:)>2))
    disp('error in region boundary')
    keyboard
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
idx=find(R>0,1);
idxList = zeros(numPixel,2);
[tx,ty]=ind2sub([dimx,dimy],idx);
idxList(1,1)=tx; idxList(1,2)=ty;

for count=2:1:numPixel
    flag=0;
    for nbr_idx=1:1:8
        px=tx+nbr(nbr_idx,1);
        py=ty+nbr(nbr_idx,2);
        if(px<1 || py<1 || px>dimx || py>dimy)
            continue;
        end
        if(R(px,py)>0)
            flag=1;
            break;
        end
    end
    
    if(flag)
        R(px,py)=0;
        idxList(count,1)=px;
        idxList(count,2)=py;
        tx=px; ty=py;
    else
        disp('error in traversing region boundary');
        keyboard;
    end
end

