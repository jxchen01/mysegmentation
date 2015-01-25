function [ctlImg,flag]=Curvature(ctlList,ctlImg,opt)

ctlImg=logical(ctlImg);
ctlLength=size(ctlList,1);


%desample into 50%
nPoints= round(ctlLength/2);
dis=[0;cumsum(sqrt(sum((ctlList(2:end,:)-ctlList(1:end-1,:)).^2,2)))];
K(:,1) = interp1(dis,ctlList(:,1),linspace(0,dis(end),nPoints));
K(:,2) = interp1(dis,ctlList(:,2),linspace(0,dis(end),nPoints));

%upsample into 200%
O(:,1)=interp(K(:,1),4);
O(:,2)=interp(K(:,2),4);

%resample back to 100%
clear K
dis=[0;cumsum(sqrt(sum((O(2:end,:)-O(1:end-1,:)).^2,2)))];
K(:,1) = interp1(dis,O(:,1),linspace(0,dis(end),ctlLength));
K(:,2) = interp1(dis,O(:,2),linspace(0,dis(end),ctlLength));

% smooth the curve
B=SnakeInternalForceMatrix2D(ctlLength,0.2,1);
S = B*K;

% compute curvature
dx=gradient(S(:,1));
ddx=gradient(dx);
dy=gradient(S(:,2));
ddy=gradient(dy);
ck=(dx.*ddy-dy.*ddx)./(hypot(dx,dy)).^3;
ck=abs(ck);

if(norm(S(1,:)-S(end,:))>2)
    buff=min([10,opt.minLength]);
    ck(1:buff)=0;
    ck(end-buff+1:1:end)=0;
end

%%%%% visual check %%%%%
% figure,imshow(ctlImg);hold on, plot(S(:,2),S(:,1),'r.');hold off

%find cutting position
cutIdx=find(ck>opt.maxCurvature);
if(numel(cutIdx)>0) 

    for i=1:1:numel(cutIdx)
        px=round(S(cutIdx(i),1));
        py=round(S(cutIdx(i),2));
        for pi=-1:1:1
            for pj=-1:1:1
                ctlImg(px+pi,py+pj)=false;
            end
        end  
    end
    flag=true;
else
    flag=false;
end
