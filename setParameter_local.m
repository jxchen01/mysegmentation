function opt=setParameter_local()

%%%% data information %%%%
sqNum=4;
numFrame=50;
idxBase=100;
imgType='.png';
fpath = 'C:\Users\jchen16\Dropbox\Private\miccai2015';

opt=struct('sqNum',sqNum,'numFrame',numFrame,...
    'idxBase',idxBase,'imgType',imgType,'filePath',fpath);

%%%% parameter for step 2 %%%%%%
opt.minArea=40; % any connected component  smaller than minArea will be removed
%opt.MBL=10;
opt.minThickness=1;
opt.maxThickness=10;

%%% paramter for step 4 %%%%
opt.maxLength=100; %%% any centerline longer than maxLength will be checked for curvature
opt.minLength=10; %%% any centerline must have more than minLength pixels
opt.maxCurvature=0.225;

%%% parameter for step 5 %%%%
opt.reconstructThickness = 4;



