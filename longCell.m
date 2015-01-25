function R=longCell(ctl, opt)


% %%%%%%%%%%%%%%%%%%%
% [Grad,~] = imgradient(Im);
% mm=max(Grad(:));
% Grad=mm-Grad;
% Grad=mat2gray(Grad)+1;
% 
% Iv= imcomplement(Im);
% Iv=mat2gray(Iv)+1;
%%%%%%%%%%%%%%%%%%%
%%% outer -medial axis %%%
%%%%%%%%%%%%%%%%%%%
R=zeros(size(ctl));
cc=bwconncomp(ctl);
labCell=labelmatrix(cc);

for i=1:1:cc.NumObjects
    ss=cc.PixelIdxList{i};
    len=length(ss);
    if(len>opt.maxLength)
        cg=ismember(labCell,i);
        sortedPixelList = SortCellPixel(cg);
        [newCtl,~]=Curvature(sortedPixelList,cg,opt);
        R(newCtl)=1;
    else
       R(ss)=1;
    end
end

% Rs=bwmorph(R,'spur');
% se_med=strel('disk',1,0);
% bg=imcomplement(imdilate(Rs,se_med));
% bg(1,:)=1;bg(end,:)=1;bg(:,1)=1;bg(:,end)=1;
% outMed = bwmorph(bg,'thin',Inf);
% 
% cc=bwconncomp(R);
% labCell=labelmatrix(cc);
% 
% voro = imcomplement(outMed);
% labVoro = labelVoronoi(labCell, voro);
% 
% R=zeros(size(ctl));
% for i=1:1:cc.NumObjects
%     ss=cc.PixelIdxList{i};
%     len=length(ss);
%     if(len>opt.maxLength)
%         rg=ismember(labVoro,i);
%         cg=ismember(labCell,i);
%         snew=optimalCut(cg, rg, Grad, Iv);
%         %snew=optimalCurve(cg, rg, Grad, Iv);
%         R(snew)=1;
%     else
%         R(ss)=1;
%     end
% end

R=bwareaopen(R, opt.minLength,8);
