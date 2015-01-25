function ctl=pruneEP(ctl)

ep=bwmorph(ctl,'endpoint');
bp=bwmorph(ctl,'branchpoint');

se=strel('square',3);
bpnbr = imdilate(bp,se);

tmp = ep & bpnbr;
ctl(tmp)=0;