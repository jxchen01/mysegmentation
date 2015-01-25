function aveThick = CellThicknesss(I,ctl)

bg=imcomplement(I);
pDist=bwdist(bg);
aveThick = mean(pDist(ctl));