function [] = dualvideo(real_data,gen_data)

structure = '%f %f %f %*s %*s';

timestamp = datestr(now,'yyyymmddTHHMMSS');
dir_ref = ['output_',timestamp];
mkdir(dir_ref);

fid = fopen(real_data);
rawdata_real = textscan(fid,structure,'Delimiter',',');
fclose(fid);

fid = fopen(gen_data);
rawdata_gen = textscan(fid,structure,'Delimiter',',');
fclose(fid);

%==Extract and Clean Data==%
data_real = cell2mat(rawdata_real);
data_real(:,1) = data_real(:,1)-data_real(1,1);
lowestID = min(min(data_real(:,2)),min(data_real(:,3)));
data_real(:,2) = data_real(:,2)-lowestID+1;
data_real(:,3) = data_real(:,3)-lowestID+1;
number_rows = size(data_real,1);
parfor i=1:number_rows
    thisrow = data_real(i,:);
    col2 = thisrow(1,2);
    col3 = thisrow(1,3);
    if col2 > col3
        thisrow(1,2) = col3;
        thisrow(1,3) = col2;
        data_real(i,:) = thisrow;
    end
end
all_IDs = [data_real(:,2); data_real(:,3)];
all_active = unique(all_IDs);
num_people = size(all_active,1);
data2 = data_real(:,2);
data3 = data_real(:,3);
for i=1:num_people
    oldID = all_active(i);
    data2(data2==oldID) = -i;
    data3(data3==oldID) = -i;
end
data_real(:,2) = -data2;
data_real(:,3) = -data3;

data_gen = cell2mat(rawdata_gen);
data_gen(:,1) = data_gen(:,1)-data_gen(1,1);
lowestID = min(min(data_gen(:,2)),min(data_gen(:,3)));
data_gen(:,2) = data_gen(:,2)-lowestID+1;
data_gen(:,3) = data_gen(:,3)-lowestID+1;
number_rows = size(data_gen,1);
parfor i=1:number_rows
    thisrow = data_gen(i,:);
    col2 = thisrow(1,2);
    col3 = thisrow(1,3);
    if col2 > col3
        thisrow(1,2) = col3;
        thisrow(1,3) = col2;
        data_gen(i,:) = thisrow;
    end
end
all_IDs = [data_gen(:,2); data_gen(:,3)];
all_active = unique(all_IDs);
num_people = size(all_active,1);
data2 = data_gen(:,2);
data3 = data_gen(:,3);
for i=1:num_people
    oldID = all_active(i);
    data2(data2==oldID) = -i;
    data3(data3==oldID) = -i;
end
data_gen(:,2) = -data2;
data_gen(:,3) = -data3;

%==Perform Analysis==%
DinT_real = deg_in_time(data_real);
DinT_gen = deg_in_time(data_gen);
create_avi(data_real,data_gen,dir_ref);
dualDDVid(DinT_real,DinT_gen,dir_ref);