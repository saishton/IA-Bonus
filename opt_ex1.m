function [ex1_opt] = opt_ex1(startex1,hrange,npoints,realFolder)

ln1 = 6.6106;
ln2 = 1.2887;

numcomp = 2;

lowex1 = max(0,startex1-hrange);
highex1 = startex1+hrange;
trials = linspace(lowex1,highex1,npoints);
comp_val = zeros(numcomp,length(trials));

RiF = ['input/',realFolder];
RtoExtract = [RiF,'/*.csv'];

RfileData = dir(RtoExtract);
RfileList = {RfileData.name};

for i=1:length(RfileList)
    currentFile = RfileList{i};
    currentClean = strrep(currentFile, '.', '');
    currentClean = strrep(currentClean, '-', '');
    currentData = pullData(RiF,currentFile,'%f %f %f %*s %*s');
    RInteractionTimes.(currentClean) = currentData.InteractionTimes_data;
end

real_fields = fieldnames(RInteractionTimes);
real_all = [];
for i = 1:numel(real_fields)
    thisdata = RInteractionTimes.(real_fields{i});
    real_all = [real_all,thisdata];
end

parfor i = 1:length(trials)
    ex1 = trials(i);
    genMatrix = model4opt(30,30000,ex1,ln1,ln2);
    gen4comp = pullData4comp(genMatrix);
    gen_data = gen4comp.InteractionTimes_data;
    comp_val(:,i) = compare4opt(real_all,gen_data); %Write me
end

[~,idx] = min(comp_val,[],2);
parfor i=1:numcomp
    ex1_opt(i) = trials(idx(i));
end