function [Ps, labmat, bw]=retrieveRegion(ctl,seg, opt)

sz=size(ctl);
se=strel('disk',opt.reconstructThickness,0);
cc=bwconncomp(ctl);
numCell=cc.NumObjects;
labmat=labelmatrix(cc);

% output variables
Ps=cell(1,numCell);
bw=zeros(sz);

for cellID=1:1:numCell
    % fetch centerline (sorted list)
    ctlPix = ismember(labmat,cellID);
    pts = SortCellPixel(ctlPix);
    
    % fetch the corresponding region

    a=bwmorph(ctlPix,'spur');
    seg_region=imdilate(a,se) & seg;
    ss=bwconncomp(seg_region,4);
    if(ss.NumObjects~=1)
        numPixels = cellfun(@numel,ss.PixelIdxList);
        [~,maxID] = max(numPixels);
        seg_region=zeros(sz);
        seg_region(ss.PixelIdxList{maxID})=1;
    else
        maxID=1;
    end
    bw(ss.PixelIdxList{maxID})=cellID;
    
    Ps{cellID}=struct('length',numel(cc.PixelIdxList{cellID}),'ctl',pts,'child',[],...
        'parent',[],'candi',[],'inflow',0,'outflow',0,'relaxinCost',0,...
        'relaxOutCost',0,'seg',seg_region);
    
end


