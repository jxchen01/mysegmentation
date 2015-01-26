opt=setParameter_local();

%%%%%%%%%%%%%%%%%%%%%%%%%
matEachFrame=cell(1,opt.numFrame);
cellEachFrame=cell(1,opt.numFrame);

for i=1:1:opt.numFrame
    disp(['frame: ',num2str(i)]);
 
    %%%%%% read raw image %%%%%%%
    raw_str=[opt.filePath,'\sq',num2str(opt.sqNum),'\raw\img0',...
        num2str(opt.idxBase+i),opt.imgType];
    I0=mat2gray(imread(raw_str));
    
    %%%% step 1: initial segmentation %%%%%%
    seg_str=[opt.filePath,'\sq',num2str(opt.sqNum),'\seg\img0',...
        num2str(opt.idxBase+i),opt.imgType];
    I1 =im2bw(imread(seg_str));
    
    %%%% step 2: remove non-cell region %%%%
    %%%% step 3: resovle intersection %%%%
    I2 = removeNonCell(I1,opt);
    
    %%%% additional check: circle %%%%%
    I2 = checkCircle(I2,opt);
    
    %%%% step 4: check long regions %%%%
    I3 = longCell(I2, opt);
    
    I4= removeCrossing(I3,opt);
    
    [cellFrame, matFrame, bw] = retrieveRegion(I4,I1,opt);
    
    %%%% save the results %%%%
    cellEachFrame{i}=cellFrame;
    matEachFrame{i}=struct('Mat',matFrame);
end

save([opt.filePath,'\sq',num2str(opt.sqNum),'\seg.mat'],'cellEachFrame','matEachFrame','-v7.3');