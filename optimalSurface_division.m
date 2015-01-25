function R=optimalSurface_division(I, inner_med, outer_med)

% input: I, grayscale image; double [0,1]
%          inner_med, binary image of inner medial axis
%          outer_med, binary image of outer medial axis
% output: R, binary segmentation

%%% parameters %%%
smoothness=2;
minJump=10;
colLength = 10;
[dimx,dimy]=size(I);

%%% class check %%%
I=mat2gray(I);
if(~islogical(inner_med))
    inner_med=(inner_med>0);
end
if(~islogical(outer_med))
    outer_med=(outer_med>0);
end

%%%% smooth %%%%
Im=normConv(I,1);

%%%% get the magnitude of gradient %%%%
[Grad,~] = imgradient(Im);
mm=max(Grad(:));
Grad=mm-Grad;
Grad=mat2gray(Grad)+1;

Iv= imcomplement(Im);
Iv=mat2gray(Iv)+1;

%%%% build the vertex-weighted graph %%%%
outerList = sortPixels(outer_med);
[xx, yy]= find(inner_med>0);
innerList = [xx, yy];
%outerList is an N-by-2 matrix. 
colWidth = size(outerList,1);
GV=zeros(colLength,colWidth);
GVtoRaw=cell(colLength,colWidth);
for i=1:1:colWidth
    tx=outerList(i,1); ty=outerList(i,2);
    idx=knnsearch(innerList,[tx,ty]);
    px=innerList(idx,1); py=innerList(idx,2);
    
    dx=(tx-px)/(colLength-1);
    dy=(ty-py)/(colLength-1);
    GV(1,i)=Grad(px,py)*Iv(px,py);
    GVtoRaw{1,i}=[px,py];
    for j=1:1:colLength-1
        sx=round(px+dx*j);
        sy=round(py+dy*j);
        GV(j+1,i)=Grad(sx,sy)*Iv(px,py);
        GVtoRaw{j+1,i}=[sx,sy];   
    end 
end


%%% convert GV to an edge-weighted graph %%%
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
%%%%% all pair shortest path %%%%%%%%
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

[dist, path, pred] = graphshortestpath(GE, ...
    numNode-1, numNode, 'Method','Acyclic');

disp(dist)

R= zeros(dimx,dimy);
len = length(path);
t0=0;
jumpFlag=0;
for k=2:1:(len-1)
    [ti,tj]=ind2sub([colLength,colWidth],path(k));
    for pp=1:1:ti
        t=GVtoRaw{pp,tj};
        R(t(1),t(2))=1;
    end
    
    if(tj>t0+1)
        disp('jump')
        disp([t0,tj])
        if(jumpFlag)
            disp('multiple jump');
            keyboard
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
        [dist_ij, path_ij, pred_ij] = graphshortestpath(GE_ij,ind1, ind2, 'Method','Acyclic');
        for kk=1:1:length(path_ij)
            [tik,tjk]=ind2sub([colLength,colWidth],path_ij(kk));
            for pp=1:1:tik
                t=GVtoRaw{pp,tjk};
                R(t(1),t(2))=1;
            end
        end
    end
    t0=tj;
end