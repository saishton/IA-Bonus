function [OnPeriods,OffPeriods] = raw(input_folder,input_filename)

minimum_to_analyse = 5;

timestamp = datestr(now,'yyyymmddTHHMMSS');
structure = '%f %f %f %*s %*s';


iF = ['input/',input_folder];
oF = ['output_',timestamp];

clean_input = strrep(input_filename, '.', '');
dir_ref = [oF,'\',clean_input];
mkdir(dir_ref);

input = [iF,'/',input_filename];

fid = fopen(input);
rawdata = textscan(fid,structure,'Delimiter',',');
fclose(fid);

%==Extract and Clean Data==%
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

%==Perform Analysis==%
DAP = node_AP(data);
EAP = edge_AP(data);

maxN = length(DAP);

%On-Off Periods:

OnPeriods = [];
OffPeriods = [];

for i=1:maxN-1
    for j=i+1:maxN
        S1 = data(:,2)==i;
        S2 = data(:,3)==j;
        T1 = data(:,3)==i;
        T2 = data(:,2)==j;
        S12 = S1 & S2;
        T12 = T1 & T2;
        ST12 = S12|T12;
        currenton = data(ST12,1);
        currentonID = (currenton./20)+1;
        onoff = zeros(1,max(data(:,1)));
        onoff(currentonID) = 1;
        switchtimes = find(diff(onoff)).*20;
        durations = diff(switchtimes);
        times1 = durations(1:2:length(durations));
        times2 = durations(2:2:length(durations));
        if onoff(1)==0
            ontimes = times1;
            offtimes = times2;
        else
            ontimes = times2;
            offtimes = times1;
        end
        
        OnPeriods = [OnPeriods,ontimes];
        OffPeriods = [OffPeriods,offtimes];
    end
end




end