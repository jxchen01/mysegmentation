function ori=OrientFit(pl,px,py)

ws = warning('off','all');  % Turn off warning

ft = fittype(['a*(x-' num2str(px) ')+',num2str(py)],...
    'dependent','y','independent','x','coefficients',{'a'});

cf = fit(pl(:,1), pl(:,2),ft);

a=confint(cf);
if(a(2)-a(1)>1)
    fty = fittype(['a*(x-' num2str(py) ')+',num2str(px)],...
    'dependent','y','independent','x','coefficients',{'a'});
    cfy = fit(pl(:,2), pl(:,1),fty);
    
    b=confint(cfy);
    if((b(2)-b(1)) < (a(2)-a(1)))
        ori=atan(1/coeffvalues(cfy));
    else
        ori=atan(coeffvalues(cf));
    end
else
    ori=atan(coeffvalues(cf));
end

warning(ws)