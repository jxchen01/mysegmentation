function labelMat = labelVoronoi(cryptLabel, voronoiRegion)

cc=bwconncomp(voronoiRegion,4);
labelMat = uint8(zeros(size(cryptLabel)));
vLabel = labelmatrix(cc);
for i=1:1:cc.NumObjects
    tmp = ismember(vLabel, i);
    tmp = (tmp>0);
    id=unique(nonzeros(uint8(tmp).*cryptLabel));
    if(numel(id)==1)
        labelMat=labelMat+uint8(id).*uint8(tmp);
    %else
    elseif(numel(id)>1)
        disp('errors');
        keyboard
    end
end
