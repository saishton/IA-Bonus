function [relativefreqs] = relfreqs(maxn,t,folder)

contact_time = 20;
format = '%f %f %f %*s %*s';

iF = ['input/',folder];
toExtract = [iF,'/*.csv'];

fileData = dir(toExtract);
fileList = {fileData.name};

freqs = zeros(1,maxn+1);

for fn=1:length(fileList)
    currentFile = fileList{fn};
    input = [iF,'/',currentFile];
    fid = fopen(input);
    rawdata = textscan(fid,format,'Delimiter',',');
    fclose(fid);
    
    data = cell2mat(rawdata);
    data(:,1) = data(:,1)-data(1,1);
    lowestID = min(min(data(:,2)),min(data(:,3)));
    data(:,2) = data(:,2)-lowestID+1;
    data(:,3) = data(:,3)-lowestID+1;
    number_rows = size(data,1);
    parfor i=1:number_rows
        thisrow = data(i,:);
        col2 = thisrow(1,2);
        col3 = thisrow(1,3);
        if col2 > col3
            thisrow(1,2) = col3;
            thisrow(1,3) = col2;
            data(i,:) = thisrow;
        end
    end
    all_IDs = [data(:,2); data(:,3)];
    all_active = unique(all_IDs);
    num_people = size(all_active,1);
    data2 = data(:,2);
    data3 = data(:,3);
    for i=1:num_people
        oldID = all_active(i);
        data2(data2==oldID) = -i;
        data3(data3==oldID) = -i;
    end
    data(:,2) = -data2;
    data(:,3) = -data3;
    
    lasttime = data(end,1);
    countermatrix = Inf(num_people);
    
    parfor i=1:num_people
        for j=1:num_people
            idx = find(data(:,2)==i & data(:,3)==j);
            if length(idx)>1
                times=data(idx,1);
                if times(1)+t<lasttime
                    times=times-times(1);
                    step = diff(times);
                    idx2 = find(step-contact_time)+1;
                    subontimes = times(idx2);
                    countermatrix(i,j) = sum(subontimes(:)<=t);
                end
            end
        end
    end
    
    thisfreqs = zeros(1,maxn+1);
    parfor i=1:maxn+1
        thisfreqs(i) = sum(countermatrix(:)==i-1);
    end
    freqs = freqs+thisfreqs;
end
relativefreqs = freqs/sum(freqs);