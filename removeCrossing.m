function R=removeCrossing(ctl, opt)

%%%% label centerline %%%%
labKernel=[1,1,1;1,0,1;1,1,1];
labMat = conv2(single(ctl), labKernel,'same');
labMat(~ctl)=0;

bp=(labMat>2);
cc=bwconncomp(bp,8);
numCross=cc.NumObjects;

if(numCross==0)
    R=ctl;
    return
end
bpIdx=labelmatrix(cc);
R=ctl;
for i=1:1:cc.NumObjects
    singleCross=ismember(bpIdx,i);
%    try
    [branchList, ~]= extractBranches(...
        singleCross,cc.PixelIdxList{i}, ctl, labMat);
%     catch
%         keyboard
%     end
    %%%% analyze the branches at current junction %%%%
    rmList = branchAnalysis(branchList, singleCross);
    R(rmList)=0;
end

%%%% final prune %%%%
R=bwmorph(R,'spur');
R=bwareaopen(R, opt.minLength, 8);
% remove the "node" caused by crossing %
R=bwmorph(R,'thin',1);





