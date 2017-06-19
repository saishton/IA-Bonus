function [opt_para] = opt_1para(high,low,npoints,realFolder)

numcomp = 2;

trials = linspace(high,low,npoints);
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
    parachange = trials(i);
    genMatrix = model4fullopt(30,30000,parachange);
    gen4comp = pullData4comp(genMatrix);
    gen_data = gen4comp.InteractionTimes_data;
    comp_val(:,i) = compare4opt(real_all,gen_data);
end

[~,idx] = min(comp_val,[],2);
parfor i=1:numcomp
    opt_para(i) = trials(idx(i));
end