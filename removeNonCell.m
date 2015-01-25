function R=removeNonCell(I, opt)

%%% remove small regions %%%
T0=bwareaopen(I,opt.minArea,4);

%%% refine the detected region %%%
T=regionRefine(T0);

%%%%%%%%%%%%%%%%%
cc=bwconncomp(T,8);
labmat=labelmatrix(cc);
R=zeros(size(I));
for i=1:1:cc.NumObjects
    reg=ismember(labmat,i);
    
    %%% extract centerline %%%
    ctl = bwmorph(reg,'thin',Inf);
    
    %%% remove extra-thin or extra-thick regions %%%
    thickness=CellThicknesss(reg,ctl);
    if(thickness<1.51 || thickness>5)
        labmat(reg)=0;
        continue;
    end
    
    %%% prune special bp and ep %%%
    % i.e. an EP in the 8-nbr of a BP 
    ctl_init = pruneEP(ctl);
    
    %%% check intersection %%%
   % try
    ctl_pruned = removeCrossing(ctl_init, opt);
%     catch
%         keyboard
%     end

    R=R | ctl_pruned;
end

R= bwareaopen(R,opt.minLength,8);