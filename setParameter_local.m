function opt=setParameter_local()

%%%% data information %%%%
sqNum=5;
numFrame=73;
idxBase=100;
imgType='.png';
fpath = 'C:\Users\jchen16\Dropbox\Private\miccai2015';

opt=struct('sqNum',sqNum,'numFrame',numFrame,...
    'idxBase',idxBase,'imgType',imgType,'filePath',fpath);

%%%% parameter for step 2 %%%%%%
opt.minArea=30;
opt.MBL=10;

%%% paramter for step 4 %%%%
opt.maxLength=70;
opt.minLength=20; %%% means len>=minLength
opt.maxCurvature=0.225;

%%% parameter for step 5 %%%%
opt.reconstructThickness = 4;



