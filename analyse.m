function [model,DataStruct] = analyse(input_folder,input_filename)

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
analysis_on = struct();
for i = 1:numel(fields_on)
    thison = OnPeriods.(fields_on{i});
    [F,X] = ecdf(thison);
    ccdf = 1-F;
    M1 = mean(thison);
    M2 = mean(thison.^2);
    M3 = mean(thison.^3);
    
    ex_lambda_start = M1;
    if ex_lambda_start<=0 || isnan(ex_lambda_start) || isinf(ex_lambda_start)
        ex_lambda_start = 0.1;
    end
    gm_b_start = (M2/M1)-M1;
    gm_a_start = M1/gm_b_start;
    if gm_b_start<=0 || isnan(gm_b_start) || isinf(gm_b_start)
        gm_b_start = 0.1;
    end
    if gm_a_start<=0 || isnan(gm_a_start) || isinf(gm_a_start)
        gm_a_start = 0.1;
    end
    rl_sigma_start = M1*sqrt(2/pi);
    if rl_sigma_start<=0 || isnan(rl_sigma_start) || isinf(rl_sigma_start)
        rl_sigma_start = 0.1;
    end
    ln_sigma_start = sqrt(log(M2*(M1^-2)));
    if ln_sigma_start<=0 || isnan(ln_sigma_start) || isinf(ln_sigma_start)
        ln_sigma_start = 0.1;
    end
    ln_mu_start = log(M1)-(0.5*ln_sigma_start^2);
    if isnan(ln_mu_start) || isinf(ln_mu_start)
        ln_mu_start = 0;
    end
    ml_beta_start = 0.5;
    ml_gamma_start = 0.5;
    [gp_k_start,gp_sigma_start,gp_theta_start] = gpSolve(M1,M2,M3);
    if isnan(gp_k_start) || isinf(gp_k_start)
        gp_k_start = 0;
    end
    if gp_sigma_start<=0 || isnan(gp_sigma_start) || isinf(gp_sigma_start)
        gp_sigma_start = 0.1;
    end
    if isnan(gp_theta_start) || isinf(gp_theta_start)
        gp_theta_start = 0;
    end
    wb_a_start = 0.5;
    wb_b_start = 0.5;
    
    try
        fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[ex_lambda_start]);
        ft_ex = fittype('expcdf(x,lambda,''upper'')','options',fo_ex);
        [cf_ex,~] = fit(X,ccdf,ft_ex);
        cv_ex = coeffvalues(cf_ex);
    catch
        cv_ex = [ex_lambda_start];
    end
    
    try
        fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[gm_a_start gm_b_start]);
        ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',fo_gm);
        [cf_gm,~] = fit(X,ccdf,ft_gm);
        cv_gm = coeffvalues(cf_gm);
    catch
        cv_gm = [gm_a_start gm_b_start];
    end
    
    try
        fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[rl_sigma_start]);
        ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',fo_rl);
        [cf_rl,~] = fit(X,ccdf,ft_rl);
        cv_rl = coeffvalues(cf_rl);
    catch
        cv_rl = [rl_sigma_start];
    end
    
    try
        fo_ln = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf 0],'Upper',[inf inf],'StartPoint',[ln_mu_start ln_sigma_start]);
        ft_ln = fittype('logncdf(x,mu,sigma,''upper'')','options',fo_ln);
        [cf_ln,~] = fit(X,ccdf,ft_ln);
        cv_ln = coeffvalues(cf_ln);
    catch
        cv_ln = [0 1];
    end
    
    try
        fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[ml_beta_start ml_gamma_start]);
        ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
        [cf_ml,~] = fit(X,ccdf,ft_ml);
        cv_ml = coeffvalues(cf_ml);
    catch
        cv_ml = [ml_beta_start ml_gamma_start];
    end
    
    try
        fo_gp = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf -inf 0],'StartPoint',[gp_k_start gp_sigma_start gp_theta_start]);
        ft_gp = fittype('gpcdf(x,k,sigma,theta,''upper'')','options',fo_gp);
        [cf_gp,~] = fit(X,ccdf,ft_gp);
        cv_gp = coeffvalues(cf_gp);
    catch
        cv_gp = [gp_k_start gp_sigma_start gp_theta_start];
    end
    
    try
        fo_wb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'StartPoint',[wb_a_start wb_b_start]);
        ft_wb = fittype('wblcdf(x,a,b,''upper'')','options',fo_wb);
        [cf_wb,~] = fit(X,ccdf,ft_wb);
        cv_wb = coeffvalues(cf_wb);
    catch
        cv_wb = [wb_a_start wb_b_start];
    end
    
    ex_lambda = cv_ex(1);
    gm_a = cv_gm(1);
    gm_b = cv_gm(2);
    rl_sigma = cv_rl(1);
    ln_mu = cv_ln(1);
    ln_sigma = cv_ln(2);
    ml_beta = cv_ml(1);
    ml_gamma = cv_ml(2);
    gp_k = cv_gp(1);
    gp_sigma = cv_gp(2);
    gp_theta = cv_gp(3);
    wb_a = cv_wb(1);
    wb_b = cv_wb(2);
    
    sortedon = sort(thison)';
    
    z_ex = expcdf(sortedon,ex_lambda);
    z_gm = gamcdf(sortedon,gm_a,gm_b);
    z_rl = raylcdf(sortedon,rl_sigma);
    z_ln = logncdf(sortedon,ln_mu,ln_sigma);
    z_ml = ones(length(sortedon),1)-mlf(ml_beta,1,-ml_gamma*sortedon.^ml_beta,6);
    z_gp = gpcdf(sortedon,gp_k,gp_sigma,gp_theta);
    z_wb = wblcdf(sortedon,wb_a,wb_b);
    
    zp_ex = exppdf(sortedon,ex_lambda);
    zp_gm = gampdf(sortedon,gm_a,gm_b);
    zp_rl = raylpdf(sortedon,rl_sigma);
    zp_ln = lognpdf(sortedon,ln_mu,ln_sigma);
    zp_ml = (-ml_beta./sortedon).*mlf(ml_beta,1,-ml_gamma*sortedon.^ml_beta,6);
    zp_gp = gppdf(sortedon,gp_k,gp_sigma,gp_theta);
    zp_wb = wblpdf(sortedon,wb_a,wb_b);
    
    stats_ex = testStatistics(sortedon,z_ex,zp_ex,20);
    stats_gm = testStatistics(sortedon,z_gm,zp_gm,20);
    stats_rl = testStatistics(sortedon,z_rl,zp_rl,20);
    stats_ln = testStatistics(sortedon,z_ln,zp_ln,20);
    stats_ml = testStatistics(sortedon,z_ml,zp_ml,20);
    stats_gp = testStatistics(sortedon,z_gp,zp_gp,20);
    stats_wb = testStatistics(sortedon,z_wb,zp_wb,20);
    
    struc_ex = struct('Scale',ex_lambda);
    struc_gm = struct('Shape',gm_a,'Scale',gm_b);
    struc_rl = struct('Scale',rl_sigma);
    struc_ln = struct('Location',ln_mu,'Scale',ln_sigma);
    struc_ml = struct('Stability',ml_beta,'Scale',ml_gamma);
    struc_gp = struct('Shape',gp_k,'Scale',gp_sigma,'Location',gp_theta);
    struc_wb = struct('Scale',wb_a,'Shape',wb_b);
    
    EX = struct('Parameters',struc_ex,'Statistics',stats_ex);
    GM = struct('Parameters',struc_gm,'Statistics',stats_gm);
    RL = struct('Parameters',struc_rl,'Statistics',stats_rl);
    LN = struct('Parameters',struc_ln,'Statistics',stats_ln);
    ML = struct('Parameters',struc_ml,'Statistics',stats_ml);
    GP = struct('Parameters',struc_gp,'Statistics',stats_gp);
    WB = struct('Parameters',struc_wb,'Statistics',stats_wb);
    
    thisStructure = struct('Exponential',EX,'Gamma',GM,'Rayleigh',RL,'LogNormal',LN,'MittagLeffler',ML,'GenPareto',GP,'Weibull',WB);
    analysis_on.(fields_on{i}) = thisStructure;
end

%Fitting for off-periods
fields_off = fieldnames(OffPeriods);
analysis_off = struct();
for i = 1:numel(fields_off)
    thisoff = OffPeriods.(fields_off{i});
    [F,X] = ecdf(thisoff);
    ccdf = 1-F;
    M1 = mean(thisoff);
    M2 = mean(thisoff.^2);
    M3 = mean(thisoff.^3);
    
    ex_lambda_start = M1;
    if ex_lambda_start<=0 || isnan(ex_lambda_start) || isinf(ex_lambda_start)
        ex_lambda_start = 0.1;
    end
    gm_b_start = (M2/M1)-M1;
    gm_a_start = M1/gm_b_start;
    if gm_b_start<=0 || isnan(gm_b_start) || isinf(gm_b_start)
        gm_b_start = 0.1;
    end
    if gm_a_start<=0 || isnan(gm_a_start) || isinf(gm_a_start)
        gm_a_start = 0.1;
    end
    rl_sigma_start = M1*sqrt(2/pi);
    if rl_sigma_start<=0 || isnan(rl_sigma_start) || isinf(rl_sigma_start)
        rl_sigma_start = 0.1;
    end
    ln_sigma_start = sqrt(log(M2*(M1^-2)));
    if ln_sigma_start<=0 || isnan(ln_sigma_start) || isinf(ln_sigma_start)
        ln_sigma_start = 0.1;
    end
    ln_mu_start = log(M1)-(0.5*ln_sigma_start^2);
    if isnan(ln_mu_start) || isinf(ln_mu_start)
        ln_mu_start = 0;
    end
    ml_beta_start = 0.5;
    ml_gamma_start = 0.5;
    [gp_k_start,gp_sigma_start,gp_theta_start] = gpSolve(M1,M2,M3);
    if isnan(gp_k_start) || isinf(gp_k_start)
        gp_k_start = 0;
    end
    if gp_sigma_start<=0 || isnan(gp_sigma_start) || isinf(gp_sigma_start)
        gp_sigma_start = 0.1;
    end
    if isnan(gp_theta_start) || isinf(gp_theta_start)
        gp_theta_start = 0;
    end
    wb_a_start = 0.5;
    wb_b_start = 0.5;
    
    try
        fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[ex_lambda_start]);
        ft_ex = fittype('expcdf(x,lambda,''upper'')','options',fo_ex);
        [cf_ex,~] = fit(X,ccdf,ft_ex);
        cv_ex = coeffvalues(cf_ex);
    catch
        cv_ex = [ex_lambda_start];
    end
    
    try
        fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[gm_a_start gm_b_start]);
        ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',fo_gm);
        [cf_gm,~] = fit(X,ccdf,ft_gm);
        cv_gm = coeffvalues(cf_gm);
    catch
        cv_gm = [gm_a_start gm_b_start];
    end
    
    try
        fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[rl_sigma_start]);
        ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',fo_rl);
        [cf_rl,~] = fit(X,ccdf,ft_rl);
        cv_rl = coeffvalues(cf_rl);
    catch
        cv_rl = [rl_sigma_start];
    end
    
    try
        fo_ln = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf 0],'Upper',[inf inf],'StartPoint',[ln_mu_start ln_sigma_start]);
        ft_ln = fittype('logncdf(x,mu,sigma,''upper'')','options',fo_ln);
        [cf_ln,~] = fit(X,ccdf,ft_ln);
        cv_ln = coeffvalues(cf_ln);
    catch
        cv_ln = [0 1];
    end
    
    try
        fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[ml_beta_start ml_gamma_start]);
        ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
        [cf_ml,~] = fit(X,ccdf,ft_ml);
        cv_ml = coeffvalues(cf_ml);
    catch
        cv_ml = [ml_beta_start ml_gamma_start];
    end
    
    try
        fo_gp = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf -inf 0],'StartPoint',[gp_k_start gp_sigma_start gp_theta_start]);
        ft_gp = fittype('gpcdf(x,k,sigma,theta,''upper'')','options',fo_gp);
        [cf_gp,~] = fit(X,ccdf,ft_gp);
        cv_gp = coeffvalues(cf_gp);
    catch
        cv_gp = [gp_k_start gp_sigma_start gp_theta_start];
    end
    
    try
        fo_wb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'StartPoint',[wb_a_start wb_b_start]);
        ft_wb = fittype('wblcdf(x,a,b,''upper'')','options',fo_wb);
        [cf_wb,~] = fit(X,ccdf,ft_wb);
        cv_wb = coeffvalues(cf_wb);
    catch
        cv_wb = [wb_a_start wb_b_start];
    end
    
    ex_lambda = cv_ex(1);
    gm_a = cv_gm(1);
    gm_b = cv_gm(2);
    rl_sigma = cv_rl(1);
    ln_mu = cv_ln(1);
    ln_sigma = cv_ln(2);
    ml_beta = cv_ml(1);
    ml_gamma = cv_ml(2);
    gp_k = cv_gp(1);
    gp_sigma = cv_gp(2);
    gp_theta = cv_gp(3);
    wb_a = cv_wb(1);
    wb_b = cv_wb(2);
    
    sortedoff = sort(thisoff)';
    
    z_ex = expcdf(sortedoff,ex_lambda);
    z_gm = gamcdf(sortedoff,gm_a,gm_b);
    z_rl = raylcdf(sortedoff,rl_sigma);
    z_ln = logncdf(sortedoff,ln_mu,ln_sigma);
    z_ml = ones(length(sortedoff),1)-mlf(ml_beta,1,-ml_gamma*sortedoff.^ml_beta,6);
    z_gp = gpcdf(sortedoff,gp_k,gp_sigma,gp_theta);
    z_wb = wblcdf(sortedoff,wb_a,wb_b);
    
    zp_ex = exppdf(sortedoff,ex_lambda);
    zp_gm = gampdf(sortedoff,gm_a,gm_b);
    zp_rl = raylpdf(sortedoff,rl_sigma);
    zp_ln = lognpdf(sortedoff,ln_mu,ln_sigma);
    zp_ml = (-ml_beta./sortedoff).*mlf(ml_beta,1,-ml_gamma*sortedoff.^ml_beta,6);
    zp_gp = gppdf(sortedoff,gp_k,gp_sigma,gp_theta);
    zp_wb = wblpdf(sortedoff,wb_a,wb_b);
    
    stats_ex = testStatistics(sortedoff,z_ex,zp_ex,20);
    stats_gm = testStatistics(sortedoff,z_gm,zp_gm,20);
    stats_rl = testStatistics(sortedoff,z_rl,zp_rl,20);
    stats_ln = testStatistics(sortedoff,z_ln,zp_ln,20);
    stats_ml = testStatistics(sortedoff,z_ml,zp_ml,20);
    stats_gp = testStatistics(sortedoff,z_gp,zp_gp,20);
    stats_wb = testStatistics(sortedoff,z_wb,zp_wb,20);
    
    struc_ex = struct('Scale',ex_lambda);
    struc_gm = struct('Shape',gm_a,'Scale',gm_b);
    struc_rl = struct('Scale',rl_sigma);
    struc_ln = struct('Location',ln_mu,'Scale',ln_sigma);
    struc_ml = struct('Stability',ml_beta,'Scale',ml_gamma);
    struc_gp = struct('Shape',gp_k,'Scale',gp_sigma,'Location',gp_theta);
    struc_wb = struct('Scale',wb_a,'Shape',wb_b);
    
    EX = struct('Parameters',struc_ex,'Statistics',stats_ex);
    GM = struct('Parameters',struc_gm,'Statistics',stats_gm);
    RL = struct('Parameters',struc_rl,'Statistics',stats_rl);
    LN = struct('Parameters',struc_ln,'Statistics',stats_ln);
    ML = struct('Parameters',struc_ml,'Statistics',stats_ml);
    GP = struct('Parameters',struc_gp,'Statistics',stats_gp);
    WB = struct('Parameters',struc_wb,'Statistics',stats_wb);
    
    thisStructure = struct('Exponential',EX,'Gamma',GM,'Rayleigh',RL,'LogNormal',LN,'MittagLeffler',ML,'GenPareto',GP,'Weibull',WB);
    analysis_off.(fields_off{i}) = thisStructure;
end

onoff2latex(analysis_on,dir_ref,'OnDistributions');
onoff2latex(analysis_off,dir_ref,'OffDistributions');

end