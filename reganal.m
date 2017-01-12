function mdl = reganal(DAP,EAP)

size = length(DAP);

ind1 = [];
ind2 = [];
dep1 = [];

for i=1:size-1
    newind1 = DAP(i+1:end);
    newind2 = ones(1,length(newind1))*DAP(i);
    newdep1 = EAP(i,i+1:end);
    ind1 = [ind1 newind1];
    ind2 = [ind2 newind2];
    dep1 = [dep1 newdep1];
end

tbl = table(ind1',ind2',dep1','VariableNames',{'AP_Vec1','AP_Vec2','AP_Edge'});

mdl = fitlm(tbl,'interactions');

end