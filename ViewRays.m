outer_length = size(outerList,1);
inner_length = size(innerList,1);
GVT=zeros(colLength,outer_length*2); % pre-assign a large size (2*length)
GVtoRaw=cell(colLength,outer_length);
colWidth = 0;
for i=1:2:outer_length
    
    tx=outerList(i,1); ty=outerList(i,2);
    
    if(i==1)
        idx=knnsearch(innerList,[tx ty]);
        px=innerList(idx,1); py=innerList(idx,2);
        %%%% determine the order of innerList to be traversed %%%
        if(idx==1)
            pre_idx=1;
        else
            ccx0=outerList(5,1);ccy0=outerList(5,2);
            cc_idx1 = mod(idx+5,inner_length);
            cc_idx2 = mod(idx-5,inner_length);
            
            ccx1=innerList(cc_idx1,1);ccy1=innerList(cc_idx1,2);
            ccx2=innerList(cc_idx2,1);ccy2=innerList(cc_idx2,2);
            
            if(hypot(ccx0-ccx1,ccy0-ccy1)> hypot(ccx0-ccx2,ccy0-ccy2))
                innerList(1:1:end,:)=innerList(end:-1:1,:);
                pre_idx = inner_length - idx +1; 
            else
                pre_idx = idx;
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        colWidth = colWidth+1;
        dx=(tx-px)/(colLength-1);
        dy=(ty-py)/(colLength-1);
        GVT(1,colWidth)=Grad(px,py)*Iv(px,py);
        GVtoRaw{1,colWidth}=[px,py];
        dbMat=zeros(dimx,dimy);
        dbMat(tx,ty)=1;
        for j=1:1:colLength-1
            sx=round(px+dx*j);
            sy=round(py+dy*j);
            GVT(j+1,colWidth)=Grad(sx,sy)*Iv(px,py);
            GVtoRaw{j+1,colWidth}=[sx,sy];
            dbMat(sx,sy)=1;
        end
        figure, imshow(dbMat)
    else
        if(pre_idx>inner_length-maxSearchAhead+1)
            sub_idx = cat(2,pre_idx:1:inner_length, 1:1:pre_idx-inner_length+maxSearchAhead-1);
        else
            sub_idx = pre_idx:1:pre_idx+maxSearchAhead-1;
        end
        
        opt_sub=knnsearch(innerList(sub_idx,:),[tx,ty]);
        opt_idx = sub_idx(opt_sub);
        
        if(opt_idx>pre_idx+4)
            for kk=pre_idx+1:2:opt_idx
                px=innerList(kk,1); py=innerList(kk,2);
                colWidth = colWidth+1;
                dx=(tx-px)/(colLength-1);
                dy=(ty-py)/(colLength-1);
                GVT(1,colWidth)=Grad(px,py)*Iv(px,py);
                GVtoRaw{1,colWidth}=[px,py];
                dbMat=zeros(dimx,dimy);
                dbMat(px,py)=1;
                for j=1:1:colLength-1
                    sx=round(px+dx*j);
                    sy=round(py+dy*j);
                    GVT(j+1,colWidth)=Grad(sx,sy)*Iv(px,py);
                    GVtoRaw{j+1,colWidth}=[sx,sy];
                    dbMat(sx,sy)=1;
                end
                figure, imshow(dbMat)
            end
        elseif((opt_idx<pre_idx && opt_idx+ inner_length > pre_idx+4 ))
            for kk=pre_idx+1:2:inner_length
                px=innerList(kk,1); py=innerList(kk,2);
                colWidth = colWidth+1;
                dx=(tx-px)/(colLength-1);
                dy=(ty-py)/(colLength-1);
                GVT(1,colWidth)=Grad(px,py)*Iv(px,py);
                GVtoRaw{1,colWidth}=[px,py];
                dbMat=zeros(dimx,dimy);
                dbMat(px,py)=1;
                for j=1:1:colLength-1
                    sx=round(px+dx*j);
                    sy=round(py+dy*j);
                    GVT(j+1,colWidth)=Grad(sx,sy)*Iv(px,py);
                    GVtoRaw{j+1,colWidth}=[sx,sy];
                    dbMat(sx,sy)=1;
                end
                dbMat(sx,sy)=1;
            end
            for kk=1:2:opt_idx
                px=innerList(kk,1); py=innerList(kk,2);
                colWidth = colWidth+1;
                dx=(tx-px)/(colLength-1);
                dy=(ty-py)/(colLength-1);
                GVT(1,colWidth)=Grad(px,py)*Iv(px,py);
                GVtoRaw{1,colWidth}=[px,py];
                dbMat=zeros(dimx,dimy);
                dbMat(px,py)=1;
                for j=1:1:colLength-1
                    sx=round(px+dx*j);
                    sy=round(py+dy*j);
                    GVT(j+1,colWidth)=Grad(sx,sy)*Iv(px,py);
                    GVtoRaw{j+1,colWidth}=[sx,sy];
                    dbMat(sx,sy)=1;
                end
                figure,imshow(dbMat)
            end
        else
            px=innerList(opt_idx,1); py=innerList(opt_idx,2);
            colWidth = colWidth+1;
            dx=(tx-px)/(colLength-1);
            dy=(ty-py)/(colLength-1);
            GVT(1,colWidth)=Grad(px,py)*Iv(px,py);
            GVtoRaw{1,colWidth}=[px,py];
            dbMat=zeros(dimx,dimy);
            dbMat(px,py)=1;
            for j=1:1:colLength-1
                sx=round(px+dx*j);
                sy=round(py+dy*j);
                GVT(j+1,colWidth)=Grad(sx,sy)*Iv(px,py);
                GVtoRaw{j+1,colWidth}=[sx,sy];
                dbMat(sx,sy)=1;
            end
            figure, imshow(dbMat)
        end
        
        pre_idx = opt_idx;
            
    end
end