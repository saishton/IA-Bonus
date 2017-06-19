function fit = testStatistics2(data1,data2,slice)

mindata = min([data1,data2]);
maxdata = max([data1,data2]);

sliced = mindata:slice:maxdata;

Phist = histcounts(data1,sliced);
Qhist = histcounts(data2,sliced);

P = Phist/sum(Phist);
Q = Qhist/sum(Qhist);

P(P==0) = 1^-50;
Q(Q==0) = 1^-50;

KL = KLDiv(P,Q);
JS = JSDiv(P,Q);

fit = struct(   'Kullback_Leibler',KL,...
                'Jensen_Shannon',JS);
end