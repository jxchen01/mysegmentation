function opt_ctlList=optimalCurve(ctl, rg, Grad, Iv)

smoothness=2;
minJump=15;
approximateLength = 2;
colLength = 10;
se_dil = strel('disk',8);
se_1 = strel('disk',1,0);
se_2 = strel('disk',2,0);
se_sq_5 = strel('square',5);
[dimx,dimy]=size(Iv);

ws = warning('off','all');  % Turn off warning

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% get outer boundary %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
dilatedCell = imdilate(ctl,se_dil);
bounding = dilatedCell & rg;
outer_region = imerode(bounding,se_1);
%outer_med = bwmorph(outer_med,'thin',1);
%outerList = sortPixels(outer_med, dimx, dimy);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% get inner boundary %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
inner_d2 = imdilate(ctl,se_2);
inner_d2(~outer_region)=false;
inner_d2=d2Refion(inner_d2);
B=bwboundaries(inner_d2,'noholes');
if(length(B)>1)
    disp('more than two objects detected');
    keyboard
end
tmp=B{1};
colWidth = size(tmp,1)-1;
inner_list(:,:) = tmp(1:end-1,:);
%inner_d1 = imerode(inner_d2,se_1);
%inner_bd = inner_d2 - inner_d1;
%inner_bd = bwmorph(inner_bd,'thin',1);
%inner_bd = pruneInnerBoundary(inner_med);
%[inner_list, colWidth] = extractInnerBoundaryList(inner_bd,dimx,dimy);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% identify region around two ends %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
new_ctl = bwmorph(inner_d2,'skel',Inf);
ep = bwmorph(new_ctl,'endpoint');
ep_region = imdilate(ep,se_sq_5);


%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% build column graph %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
GV=zeros(colLength,colWidth);
GVtoRaw=cell(colLength,colWidth);
botV=ones(1,colWidth)*colLength;
for i=1:1:colWidth
    tx=inner_list(i,1);
    ty=inner_list(i,2);
    GVtoRaw{1,i}=[tx,ty];
    GV(1,i)=Grad(tx,ty)*Iv(tx,ty);
    
    %%%% get the normal direction %%%%%
    if(ep_region(tx,ty))        
        if(i==1)
            ddx=(inner_list(2,1)-inner_list(end,1))/2;
            ddy=(inner_list(2,2)-inner_list(end,2))/2;
        elseif(i==colWidth)
            ddx=(inner_list(1,1)-inner_list(end-1,1))/2;
            ddy=(inner_list(1,2)-inner_list(end-1,2))/2;
        else
            ddx=(inner_list(i+1,1)-inner_list(i-1,1))/2;
            ddy=(inner_list(i+1,2)-inner_list(i-1,2))/2;
        end
        nn=hypot(ddx,ddy);
        dx=-ddy/nn; dy=ddx/nn;
        
    else   
        if(i<=approximateLength)
            pl=cat(1,inner_list(end-approximateLength+i:1:end,:),...
                inner_list(1:1:i+approximateLength,:));
        elseif(i>colWidth-approximateLength)
            pl=cat(1,inner_list(i-approximateLength:1:end,:),...
                inner_list(1:1:approximateLength-colWidth+i,:));
        else
            pl=inner_list(i-approximateLength:1:i+approximateLength,:);
        end
        
        ori = OrientFit(pl,tx,ty);
        if(abs(ori-pi/2)<0.001)
            dx=1;dy=0;
        else
            dk=tan(ori); nn=sqrt(1+dk*dk);
            dx=dk/nn;dy=-1/nn;
        end
    end
    %%%%%% direction of the ray %%%%%%
    shootDirection=0;
    for j=1:1:colLength
        px=round(tx+dx*j);
        py=round(ty+dy*j);
        if(px~=tx || py~=ty)
            if(px>0 && px<=dimx && py>0 && py<=dimy && inner_d2(px,py))
                shootDirection=-1;
                break;
            end
        end
        
        px=round(tx-dx*j);
        py=round(ty-dy*j);
        if(px~=tx||py~=ty)
            if(px>0 && px<=dimx && py>0 && py<=dimy && inner_d2(px,py))
                shootDirection=1;
                break;
            end
        end
    end
    if(shootDirection==0)
        disp('error in shooting the ray');
        keyboard;
    end
    
    %%%%%%%%% shoot a ray %%%%%%%%%%%%%
    for j=1:1:colLength-1
        tj=shootDirection*j;
        px=round(tx+tj*dx);
        py=round(ty+tj*dy);
        if(px<1 || px>dimx || py<1 || py>dimy)
            continue;
        end
        if(outer_region(px,py))
            GVtoRaw{j+1,i}=[px,py];
            GV(j+1,i)=Grad(px,py)*Iv(px,py);
        else
            botV(i)=j+1;
            break;
        end
    end
    
end

%%% convert G to an edge-weighted graph %%%
numNode=colWidth*colLength+2;
GE_init=zeros(numNode*(2*smoothness+1),3);
sid=0;
for i=1:1:(colWidth-1)
    for j=1:1:botV(i)
        ind1=sub2ind([colLength,colWidth],j,i);
        for k=-smoothness:1:smoothness
            j2=j+k;
            if(j2<1 || j2>botV(i+1))
                continue;
            end
            ind2=sub2ind([colLength,colWidth],j2,i+1);
            sid=sid+1;
            GE_init(sid,1)=ind1;
            GE_init(sid,2)=ind2;
            GE_init(sid,3)=GV(ind2);
        end
    end
end

sid0=sid;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% all pair shortest path %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

GE_all= sparse(GE_init(1:sid,1), GE_init(1:sid,2), GE_init(1:sid,3), numNode-2, numNode-2);
dist_all = graphallshortestpaths(GE_all);
rangeMat = zeros(colWidth,colWidth);
for i=1:1:colWidth
    ind1=sub2ind([colLength,colWidth],1,i);
    for j=(i+minJump):1:colWidth
        ind2=sub2ind([colLength,colWidth],1,j);
        rangeMat(i,j)=dist_all(ind1,ind2)/(j-i);
    end
end
rangeMat(rangeMat==0)=Inf;

for i=2:1:colWidth
    ind1=sub2ind([colLength,colWidth],1,i);
    for j=(i+minJump+2):1:colWidth
        ind2=sub2ind([colLength,colWidth],1,j);
        %queryMat=rangeMat((i+1):1:(j-2),(i+1):1:(j-1));
        %minCost = min(queryMat(:));
        %if(minCost<Inf)
        if(rangeMat(i+1,j-1)<Inf)
            sid=sid+1;
            GE_init(sid,1)=ind1;
            GE_init(sid,2)=ind2;
            %GE_init(sid,3)=(j-i+1)*minCost;
            GE_init(sid,3)=(j-i+1)*rangeMat(i+2,j-2);
        end
    end
end

%%%% node S %%%%
for j=1:1:botV(1)
    ind2=sub2ind([colLength,colWidth],j,1);
    sid=sid+1;
    GE_init(sid,1)=numNode-1;
    GE_init(sid,2)=ind2;
    GE_init(sid,3)=GV(ind2);
end

%%%% node T %%%%
for j=1:1:botV(colWidth)
    ind1=sub2ind([colLength,colWidth],j,colWidth);
    sid=sid+1;
    GE_init(sid,1)=ind1;
    GE_init(sid,2)=numNode;
    GE_init(sid,3)=1;
end

%%%% the final graph %%%%
GE= sparse(GE_init(1:sid,1), GE_init(1:sid,2), GE_init(1:sid,3), numNode, numNode);

[~, path, ~] = graphshortestpath(GE, ...
    numNode-1, numNode, 'Method','Acyclic');

%disp(dist)

R= zeros(dimx,dimy);
len = length(path);
jumpFlag=0;

[ti,tj]=ind2sub([colLength,colWidth],path(2));
t=GVtoRaw{ti,tj};
if(~isempty(t))
    R(t(1),t(2))=1;
else
    disp('wrong path');
    keyboard
end
px=t(1);py=t(2);
px0=t(1);py0=t(2);
t0=tj;

for k=3:1:(len-1)
    [ti,tj]=ind2sub([colLength,colWidth],path(k));
    t=GVtoRaw{ti,tj};
    [xx,yy]=bresenham(px,py,t(1),t(2));
    
    for j=1:1:length(xx)
        R(xx(j), yy(j))=1;
    end
    px=t(1);py=t(2);

    if(tj>t0+1)
        %disp('jump')
        %disp([t0,tj])
        if(jumpFlag)
            %disp('multiple jump');
            %keyboard
            opt_ctlList = find(ctl>0);
            return
        end
        jumpFlag=1;
        %queryMat=rangeMat((t0+1):1:(t0+10),(tj-10):1:(tj-1));
        %queryMat=rangeMat((t0+1):1:(tj-1),(t0+1):1:(tj-1));
        %[tmp, a] = min(queryMat);
        %[tmp, b] = min(tmp);
        %disp([t0+a(b),t0+b]);
        ind1=sub2ind([colLength,colWidth],1, t0+2);
        ind2=sub2ind([colLength,colWidth],1, tj-2);
        GE_ij= sparse(GE_init(1:sid0,1), GE_init(1:sid0,2), GE_init(1:sid0,3), numNode-2, numNode-2);
        [~, path_ij, ~] = graphshortestpath(GE_ij,ind1, ind2, 'Method','Acyclic');
        
        [tik,tjk]=ind2sub([colLength,colWidth],path_ij(1));
        t=GVtoRaw{tik,tjk};
        pxk=t(1); pyk=t(2);
        for kk=2:1:length(path_ij)
            [tik,tjk]=ind2sub([colLength,colWidth],path_ij(kk));
            t=GVtoRaw{tik,tjk};
            [xx,yy]=bresenham(pxk,pyk,t(1),t(2)); 
            for j=1:1:length(xx)
                R(xx(j), yy(j))=1;
            end
            pxk=t(1);pyk=t(2);
        end
    end
    
    t0=tj;
end
%%% link the first and the last node %%%
[xx,yy]=bresenham(px0,py0,px,py);
for j=1:1:length(xx)
    R(xx(j), yy(j))=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% extract the centerline %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
Rs=imfill(R,'holes');
opt_ctl=bwmorph(Rs,'thin',Inf);
opt_ctlList = find(opt_ctl>0);

warning(ws)

end

function R=d2Refion(R)
%%% remove small spur/bump on boundary %%%
nbrKernel=[0,1,0;1,0,1;0,1,0];
resp=conv2(single(R),single(nbrKernel),'same');
resp(~R)=0;
R(resp<2)=0;
end

function [pList, numPixel]=extractInnerBoundaryList(innerBoundary,dimx,dimy)
    nbr=[1,1;1,0;1,-1;0,1;0,-1;-1,1;-1,0;-1,-1];
    
    innerBoundary=int8(innerBoundary);
    idx=find(innerBoundary>0);
    [tx,ty]=ind2sub([dimx,dimy],idx(1));
    
    numPixel = numel(idx);
    pList=zeros(numPixel,2);
    sig=1;
    
    
    pList(sig,1)=tx; pList(sig,2)=ty;
    innerBoundary(tx,ty)=-1;
    
    flag=1;
    while(flag)
        flag=0;
        for i=1:1:8
            px=tx+nbr(i,1);
            py=ty+nbr(i,2);
            if(px<1 || px>dimx || py<1 || py>dimy)
                continue;
            end
            if(innerBoundary(px,py)>0)
                sig=sig+1;
                pList(sig,1)=px;pList(sig,2)=py;
                innerBoundary(px,py)=-1;
                tx=px;ty=py;
                flag=1;
                break;
            end
        end
    end
    if(sig~=numPixel)
        disp('error in extracting inner boundary');
        keyboard;
    end
end