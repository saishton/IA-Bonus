function [ex1,ln1,ln2] = analyse3(input_folder,input_filename)

minimum_to_analyse = 1;

structure = '%f %f %f %*s %*s';


iF = ['input/',input_folder];
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
%Activity Potentials and Degree-in-Time:
DAP = node_AP(data);
EAP = edge_AP(data);
%model = reganal(DAP,EAP);
%DinT = deg_in_time(data);
%changes = deg_change(DinT);
%DDVid(DinT,dir_ref);

maxN = length(DAP);

% %Degree-in-Time Changes
% DataStruct = struct();
%
% for i=1:maxN
%     changeHist = figure();
%     histogram(changes(i,:),length(unique(changes(i,:))))
%     title_str = sprintf('DEGREE CHANGE HISTOGRAM\nNode-AP: %d', DAP(i));
%     title(title_str);
%     imagefilename = [dir_ref,'/node-',num2str(i),'_degree-change-hist.png'];
%     print(imagefilename,'-dpng')
%     close(changeHist);
%
%     x = [0:20:max(data(:,1))];
%     y = DinT(i,:);
%     DinTfig = figure();
%     plot(x,y);
%     title_str = sprintf('DEGREE CHANGE IN TIME\nNode-AP: %d', DAP(i));
%     title(title_str);
%     imagefilename = [dir_ref,'/node-',num2str(i),'_degree-change-in-time.png'];
%     print(imagefilename,'-dpng')
%     close(DinTfig);
%
%     [F,X] = ecdf(changes(i,:));
%     ccdf = 1-F;
%     M1 = mean(changes(i,:));
%     M2 = mean(changes(i,:).^2);
%     mu_g = M1;
%     mu_s = sqrt(M2-M1^2);
%     fo = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf 0],'Upper',[inf inf],'StartPoint',[mu_g mu_s]);
%     ft = fittype('normcdf(x,mu,sigma,''upper'')','options',fo);
%     [cf,~] = fit(X,ccdf,ft);
%     cv = coeffvalues(cf);
%     n_mu = cv(1);
%     n_si = cv(2);
%     ccdf_n = normcdf(X,n_mu,n_si,'upper');
%
%     fitting = figure();
%     hold on
%     plot(X,ccdf,'o');
%     plot(X,ccdf_n);
%     title_str = sprintf('DEGREE CHANGE CCDF\nNode-AP: %d', DAP(i));
%     title(title_str);
%     hold off
%     imagefilename = [dir_ref,'/node-',num2str(i),'_degree-change-ccdf.png'];
%     print(imagefilename,'-dpng')
%     close(fitting);
%
%     thisparas = struct('EAP',EAP(i),'mu',n_mu,'sigma_squared',n_si^2);
%
%     z = normcdf(sort(changes(i,:))',n_mu,n_si);
%     zp = normpdf(sort(changes(i,:))',n_mu,n_si);
%     thisstats = testStatistics(sort(changes(i,:))',z,zp,0);
%
%     paras_ref = sprintf('PARAS_N%d', i);
%     stats_ref = sprintf('STATS_N%d', i);
%
%     DataStruct.(paras_ref) = thisparas;
%     DataStruct.(stats_ref) = thisstats;
% end

%On-Off Periods:

OnPeriods = struct();
OffPeriods = struct();

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
        
        ID_ref = sprintf('n%d_n%d', i,j);
        if length(ontimes)>=minimum_to_analyse && length(unique(ontimes))>=3
            OnPeriods.(ID_ref) = ontimes;
        end
        if length(offtimes)>=minimum_to_analyse && length(unique(ontimes))>=3
            OffPeriods.(ID_ref) = offtimes;
        end
    end
end

%Fitting for on-periods
fields_on = fieldnames(OnPeriods);
ex1 = [];
for i = 1:numel(fields_on)
    thison = OnPeriods.(fields_on{i});
    M1 = mean(thison);
    ex1 = [ex1,M1];
end

ln1 = 0;
ln2 = 0;
end