function opt_ctlList=optimalCut(ctl, rg, Grad, Iv)

smoothness=2;
minJump=10;
colLength = 10;
maxSearchAhead = 20;
[dimx,dimy]=size(Iv);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% get outer boundary %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
se_dil = strel('disk',8);
%se_1 = strel('disk',1);
dilatedCell = imdilate(ctl,se_dil);
bounding = dilatedCell & rg;

%outer_med = imdilate(bounding, se_1) - bounding;
%outer_med = bwmorph(outer_med,'thin',1);
%outer_ep = bwmorph(outer_med,'endpoint');

B = bwboundaries(bounding,'noholes');
if(numel(B)~=1)
    disp('error in identifying bounderies');
    keyboard;
end
tmp=B{1};
outerList(:,:)=tmp(1:end-1, :);
clear tmp B
%outerList = sortPixels(outer_med, dimx, dimy);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% get inner boundary %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
ctl=bwmorph(ctl,'spur');
ctl=bwmorph(ctl,'spur');

BC=bwboundaries(ctl,'noholes');
tmp=BC{1};
innerList(:,:) = tmp(1:end-1,:);
%innerList = SortCellPixel(ctl);
%[inX,inY] = ind2sub([dimx,dimy], ctlList(2:1:(end-1)));
%innerList=cat(2,inX, inY);
clear tmp BC

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% build column graph %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
outer_length = size(outerList,1);
inner_length = size(innerList,1);
GVT=zeros(colLength,outer_length*2); % pre-assign a large size (2*length)
GVtoRaw=cell(colLength,outer_length);
colWidth = 0;
for i=1:2:outer_length
    tx=outerList(i,1); ty=outerList(i,2);
    if(i==1)
        idx=knnsearch(innerList,[tx ty]);
        px=innerList(idx,1); py=innerList(idx,2);
        
        %%%% determine the order of innerList to be traversed %%%
        if(idx==1)
            pre_idx=1;
        else
            ccx0=outerList(5,1);ccy0=outerList(5,2);
            cc_idx1 = mod(idx+5,inner_length);
            if(cc_idx1==0)
                cc_idx1=inner_length;
            end
            cc_idx2 = mod(idx-5,inner_length);
            if(cc_idx2==0)
                cc_idx2=inner_length;
            end
            
            ccx1=innerList(cc_idx1,1);ccy1=innerList(cc_idx1,2);
            ccx2=innerList(cc_idx2,1);ccy2=innerList(cc_idx2,2);
            
            if(hypot(ccx0-ccx1,ccy0-ccy1)> hypot(ccx0-ccx2,ccy0-ccy2))
                innerList(1:1:end,:)=innerList(end:-1:1,:);
                pre_idx = inner_length - idx +1; 
            else
                pre_idx = idx;
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        colWidth = colWidth+1;
        dx=(tx-px)/(colLength-1);
        dy=(ty-py)/(colLength-1);
        GVT(1,colWidth)=Grad(px,py)*Iv(px,py);
        GVtoRaw{1,colWidth}=[px,py];
        for j=1:1:colLength-1
            sx=round(px+dx*j);
            sy=round(py+dy*j);
            GVT(j+1,colWidth)=Grad(sx,sy)*Iv(px,py);
            GVtoRaw{j+1,colWidth}=[sx,sy];
        end
    else
        if(pre_idx>inner_length-maxSearchAhead+1)
            sub_idx = cat(2,pre_idx:1:inner_length, 1:1:pre_idx-inner_length+maxSearchAhead-1);
        else
            sub_idx = pre_idx:1:pre_idx+maxSearchAhead-1;
        end
        
        opt_sub=knnsearch(innerList(sub_idx,:),[tx,ty]);
        opt_idx = sub_idx(opt_sub);
        
        if(opt_idx>pre_idx+4)
            for kk=pre_idx+1:2:opt_idx
                px=innerList(kk,1); py=innerList(kk,2);
                colWidth = colWidth+1;
                dx=(tx-px)/(colLength-1);
                dy=(ty-py)/(colLength-1);
                GVT(1,colWidth)=Grad(px,py)*Iv(px,py);
                GVtoRaw{1,colWidth}=[px,py];
                for j=1:1:colLength-1
                    sx=round(px+dx*j);
                    sy=round(py+dy*j);
                    GVT(j+1,colWidth)=Grad(sx,sy)*Iv(px,py);
                    GVtoRaw{j+1,colWidth}=[sx,sy];
                end
            end
        elseif((opt_idx<pre_idx && opt_idx+ inner_length > pre_idx+4 ))
            for kk=pre_idx+1:2:inner_length
                px=innerList(kk,1); py=innerList(kk,2);
                colWidth = colWidth+1;
                dx=(tx-px)/(colLength-1);
                dy=(ty-py)/(colLength-1);
                GVT(1,colWidth)=Grad(px,py)*Iv(px,py);
                GVtoRaw{1,colWidth}=[px,py];
                for j=1:1:colLength-1
                    sx=round(px+dx*j);
                    sy=round(py+dy*j);
                    GVT(j+1,colWidth)=Grad(sx,sy)*Iv(px,py);
                    GVtoRaw{j+1,colWidth}=[sx,sy];
                end
            end
            for kk=1:2:opt_idx
                px=innerList(kk,1); py=innerList(kk,2);
                colWidth = colWidth+1;
                dx=(tx-px)/(colLength-1);
                dy=(ty-py)/(colLength-1);
                GVT(1,colWidth)=Grad(px,py)*Iv(px,py);
                GVtoRaw{1,colWidth}=[px,py];
                for j=1:1:colLength-1
                    sx=round(px+dx*j);
                    sy=round(py+dy*j);
                    GVT(j+1,colWidth)=Grad(sx,sy)*Iv(px,py);
                    GVtoRaw{j+1,colWidth}=[sx,sy];
                end
            end
        else
            px=innerList(opt_idx,1); py=innerList(opt_idx,2);
            colWidth = colWidth+1;
            dx=(tx-px)/(colLength-1);
            dy=(ty-py)/(colLength-1);
            GVT(1,colWidth)=Grad(px,py)*Iv(px,py);
            GVtoRaw{1,colWidth}=[px,py];
            for j=1:1:colLength-1
                sx=round(px+dx*j);
                sy=round(py+dy*j);
                GVT(j+1,colWidth)=Grad(sx,sy)*Iv(px,py);
                GVtoRaw{j+1,colWidth}=[sx,sy];
            end
        end
        
        pre_idx = opt_idx;
            
    end
end

GV(:,:)=GVT(:,1:1:colWidth);

% colWidth = size(outerList,1);
% GV=zeros(colLength,colWidth);
% GVtoRaw=cell(colLength,colWidth);
% for i=1:1:colWidth
%     tx=outerList(i,1); ty=outerList(i,2);
%     if(i==1)
%         idx=knnsearch(innerList,[tx ty]);
%         px=innerList(idx,1); py=innerList(idx,2);
%         pre_idx = idx;
%     else
%         if(pre_idx>colWidth-maxSearchAhead+1)
%             sub_idx = cat(2,pre_idx:1:colWidth, 1:1:pre_idx-colWidth+maxSearchAhead-1);
%         else
%             sub_idx = pre_idx:1:pre_idx+maxSearchAhead-1;
%         end
%         opt_sub=knnsearch(innerList(sub_idx,:),[tx,ty]);
%         opt_idx = sub_idx(opt_sub);
%         
%     end
%     
%     dx=(tx-px)/(colLength-1);
%     dy=(ty-py)/(colLength-1);
%     GV(1,i)=Grad(px,py)*Iv(px,py);
%     GVtoRaw{1,i}=[px,py];
%     for j=1:1:colLength-1
%         sx=round(px+dx*j);
%         sy=round(py+dy*j);
%         GV(j+1,i)=Grad(sx,sy)*Iv(px,py);
%         GVtoRaw{j+1,i}=[sx,sy];
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% convert GV to an edge-weighted graph %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

numNode=colWidth*colLength+2;
GE_init=zeros(numNode*(2*smoothness+1),3);
sid=0;
for i=1:1:(colWidth-1)
    for j=1:1:colLength
        ind1=sub2ind([colLength,colWidth],j,i);
        for k=-smoothness:1:smoothness
            j2=j+k;
            if(j2<1 || j2>colLength)
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
for j=1:1:colLength
    ind2=sub2ind([colLength,colWidth],j,1);
    sid=sid+1;
    GE_init(sid,1)=numNode-1;
    GE_init(sid,2)=ind2;
    GE_init(sid,3)=GV(ind2);
end

%%%% node T %%%%
for j=1:1:colLength
    ind1=sub2ind([colLength,colWidth],j,colWidth);
    sid=sid+1;
    GE_init(sid,1)=ind1;
    GE_init(sid,2)=numNode;
    GE_init(sid,3)=1;
end

% sid=sid+1;
% GE_init(sid,1)=numNode;
% GE_init(sid,2)=numNode-1;
% GE_init(sid,3)=1;

GE= sparse(GE_init(:,1), GE_init(:,2), GE_init(:,3), numNode, numNode);

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
        disp('jump')
        disp([t0,tj])
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
if(jumpFlag>0)
    cg=imfill(R,'holes');
    Rs=imcomplement(R);
    Rs(~cg)=0;
else
    Rs=imfill(R,'holes');
end

idx=find(R>0);
[tx,ty]=ind2sub([dimx,dimy],idx);
for i=1:1:length(idx)
    ss= Rs(tx(i),ty(i)+1)+ Rs(tx(i),ty(i)-1) + Rs(tx(i)+1,ty(i)) + Rs(tx(i)-1,ty(i));
    if(ss<2)
        Rs(tx(i),ty(i))=0;
    end
end

opt_ctl=bwmorph(Rs,'thin',Inf);
opt_ctlList = find(opt_ctl>0);

