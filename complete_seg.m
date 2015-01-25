
opt=setParameter_local();

%%%%%%%%%%%%%%%%%%%%%%%%%
img=cell(1,opt.numFrame);
for i=30:1:opt.numFrame
    disp(['frame: ',num2str(i)]);
    str=['../sq',num2str(opt.sqNum),'/raw/img0',...
        num2str(opt.idxBase+i),opt.imgType];
 
    %%%%%% read raw image %%%%%%%
    I0=imread(str);
    
    %%%% step 1: initial segmentation %%%%%%
    I1 = prediction(myModel, I0, numPixel, opt);
    
    %%%% step 2: remove non-cell region %%%%
    %%%% step 3: resovle intersection %%%%
    I2 = removeNonCell(I1,opt);
    
    %%%% additional check: circle %%%%%
    I2 = checkCircle(I2,opt);
    
    %%%% step 4: check long regions %%%%
    I3 = longCell(I2, opt);
    
    I4= removeCrossing(I3,opt);
    
    [bw,Ps] = retrieveRegion(I4,I1,opt);
    
    %%%% save the results %%%%
    img{i}=struct('Raw',I0,'preseg',I1,'refined',I2,'finalBW',I4);
end

save(['../sq',num2str(opt.sqNum),'/data.mat'],'img','-v7.3');