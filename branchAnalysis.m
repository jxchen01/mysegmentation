function rmList = branchAnalysis(branchList, singleCross)

[dimx,dimy]=size(singleCross);
numBranch=numel(branchList);

%%% more than three branches %%%
if(numBranch>3)
    minLinear=pi;
    b1=0;b2=0;
    for i=1:1:numBranch-1
        for j=(i+1):1:numBranch
            ta=branchList{i}.orientation;
            tb=branchList{j}.orientation;
            tt=abs(ta-tb);
            tmp=min([tt, abs( tt- pi)]);
            if(tmp<minLinear)
                minLinear=tmp;
                b1=i;b2=j;
            end
        end
    end
    rmList=[];
    if(minLinear<0.52) % co-linear exists, pi/9=0.52
        rmIdx=setdiff([1:1:numBranch],[b1,b2]);
        for k=1:1:numel(rmIdx)
            br=branchList{rmIdx(k)}.pixelList;
        	len=branchList{rmIdx(k)}.branchLength;
            if(len>1)
                rmList=cat(1,rmList,  sub2ind([dimx,dimy],br(1:2,1), br(1:2,2)));
            else
                rmList=cat(1, rmList, sub2ind([dimx,dimy],br(1,1), br(1,2)));
            end
        end
    else % no co-linear
        rmList=find(singleCross>0);
    end

else
    %%%% exactly three branches %%%%
    %%% type index:                     %%%%
    %%% 1. major branches,          %%%%
    %%% 2. minor branches,
    %%% 3. connection branches   %%%%
    %%%%%%%%%%%%%%%%%%%%
    idxList=[1,2;1,3;2,3];
    
    % case 1: co-linear pair exists,
    %            cut the remaining one with 2 pixels
    % case 2: no co-linear pair,
    %            remove the junction pixels
    
    minLinear=pi;
    minIdx=0;
    for i=1:1:3
        ta=branchList{idxList(i,1)}.orientation;
        tb=branchList{idxList(i,2)}.orientation;
        tt=abs(ta-tb);
        tmp=min([tt, abs( tt- pi)]);
        if(tmp<minLinear)
            minLinear=tmp;
            minIdx=i;
        end
    end
    
    if(minLinear<0.52) % co-linear exists, pi/9=0.52
        rmIdx=setdiff([1,2,3],idxList(minIdx,:));
        br=branchList{rmIdx}.pixelList;
        len=branchList{rmIdx}.branchLength;
        if(len>1)
            rmList=sub2ind([dimx,dimy],br(1:2,1), br(1:2,2));
        else
            rmList=sub2ind([dimx,dimy],br(1,1), br(1,2));
        end
    else % no co-linear
        rmList=find(singleCross>0);
    end


end


