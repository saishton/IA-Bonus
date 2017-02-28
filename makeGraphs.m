function [] = makeGraphs(gen_data,real_data,variable,dir_ref)

cleanName = strrep(variable, '.', '');
cleanName = strrep(cleanName, '-', '');

gen_fields = fieldnames(gen_data);
real_fields = fieldnames(real_data);

plotallthethings = figure();
hold on
for i = 1:numel(gen_fields)
    thisdata = gen_data.(gen_fields{i});
    [F,X] = ecdf(thisdata);
    ccdf = 1-F;
    plot(X,ccdf,'o');
end
for i = 1:numel(real_fields)
    thisdata = real_data.(real_fields{i});
    [F,X] = ecdf(thisdata);
    ccdf = 1-F;
    plot(X,ccdf,'x');
end
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(variable);
ylabel('CCDF');
imagefilename1 = [dir_ref,'/',cleanName,'_indiv.png'];
print(imagefilename1,'-dpng')
hold off
close(plotallthethings);

gen_all = [];
real_all = [];
for i = 1:numel(gen_fields)
    thisdata = gen_data.(gen_fields{i});
    gen_all = [gen_all,thisdata];
end
for i = 1:numel(real_fields)
    thisdata = real_data.(real_fields{i});
    real_all = [real_all,thisdata];
end

plotsomeofthethings = figure();
hold on
[Fg,Xg] = ecdf(gen_all);
ccdfg = 1-Fg;
plot(Xg,ccdfg,'o');
[Fr,Xr] = ecdf(real_all);
ccdfr = 1-Fr;
plot(Xr,ccdfr,'x');
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel(variable);
ylabel('CCDF');
imagefilename2 = [dir_ref,'/',cleanName,'_combine.png'];
print(imagefilename2,'-dpng')
hold off
close(plotsomeofthethings);
end