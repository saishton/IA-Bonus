function [] = prepModel(input_folder,ondist,offdist)

minimum_to_analyse = 5;

timestamp = datestr(now,'yyyymmddTHHMMSS');
structure = '%f %f %f %*s %*s';

iF = ['input/',input_folder];
oF = ['output_',timestamp];

dir_ref = [oF,'\',input_folder];
mkdir(dir_ref);

toExtract = [iF,'/*.csv'];
fileData = dir(toExtract);
fileList = {fileData.name};

if strcmp(ondist,'Exponential')
    paracount_on = 1;
elseif strcmp(ondist,'Gamma')
    paracount_on = 2;
elseif strcmp(ondist,'Rayleigh')
    paracount_on = 1;
elseif strcmp(ondist,'LogNormal')
    paracount_on = 2;
elseif strcmp(ondist,'MittagLeffler')
    paracount_on = 2;
elseif strcmp(ondist,'GenPareto')
    paracount_on = 3;
elseif strcmp(ondist,'Weibull')
    paracount_on = 2;
end
if strcmp(offdist,'Exponential')
    paracount_off = 1;
elseif strcmp(offdist,'Gamma')
    paracount_off = 2;
elseif strcmp(offdist,'Rayleigh')
    paracount_off = 1;
elseif strcmp(offdist,'LogNormal')
    paracount_off = 2;
elseif strcmp(offdist,'MittagLeffler')
    paracount_off = 2;
elseif strcmp(offdist,'GenPareto')
    paracount_off = 3;
elseif strcmp(offdist,'Weibull')
    paracount_off = 2;
end

para_on = [];
para_off = [];

for file=1:length(fileList)
    currentFile = fileList{file};
    input = [iF,'/',currentFile];
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
    
    %On-Off Periods:
    DAP = node_AP(data);
    maxN = length(DAP);
    
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
    
    thispara_on = zeros(numel(fields_on),paracount_on);
    
    
    for i = 1:numel(fields_on)
        thison = OnPeriods.(fields_on{i});
        [F,X] = ecdf(thison);
        ccdf = 1-F;
        M1 = mean(thison);
        M2 = mean(thison.^2);
        M3 = mean(thison.^3);
        
        if strcmp(ondist,'Exponential')
            ex_lambda_start = M1;
            if ex_lambda_start<=0 || isnan(ex_lambda_start) || isinf(ex_lambda_start)
                ex_lambda_start = 0.1;
            end
            try
                fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[ex_lambda_start]);
                ft_ex = fittype('expcdf(x,lambda,''upper'')','options',fo_ex);
                [cf_ex,~] = fit(X,ccdf,ft_ex);
                thispara_on(i,:) = coeffvalues(cf_ex);
            catch
                thispara_on(i,:) = [ex_lambda_start];
            end
        elseif strcmp(ondist,'Gamma')
            gm_b_start = (M2/M1)-M1;
            gm_a_start = M1/gm_b_start;
            if gm_b_start<=0 || isnan(gm_b_start) || isinf(gm_b_start)
                gm_b_start = 0.1;
            end
            if gm_a_start<=0 || isnan(gm_a_start) || isinf(gm_a_start)
                gm_a_start = 0.1;
            end
            try
                fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[gm_a_start gm_b_start]);
                ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',fo_gm);
                [cf_gm,~] = fit(X,ccdf,ft_gm);
                thispara_on(i,:) = coeffvalues(cf_gm);
            catch
                thispara_on(i,:) = [gm_a_start gm_b_start];
            end
        elseif strcmp(ondist,'Rayleigh')
            rl_sigma_start = M1*sqrt(2/pi);
            if rl_sigma_start<=0 || isnan(rl_sigma_start) || isinf(rl_sigma_start)
                rl_sigma_start = 0.1;
            end
            try
                fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[rl_sigma_start]);
                ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',fo_rl);
                [cf_rl,~] = fit(X,ccdf,ft_rl);
                thispara_on(i,:) = coeffvalues(cf_rl);
            catch
                thispara_on(i,:) = [rl_sigma_start];
            end
        elseif strcmp(ondist,'LogNormal')
            ln_sigma_start = sqrt(log(M2*(M1^-2)));
            if ln_sigma_start<=0 || isnan(ln_sigma_start) || isinf(ln_sigma_start)
                ln_sigma_start = 0.1;
            end
            ln_mu_start = log(M1)-(0.5*ln_sigma_start^2);
            if isnan(ln_mu_start) || isinf(ln_mu_start)
                ln_mu_start = 0;
            end
            try
                fo_ln = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf 0],'Upper',[inf inf],'StartPoint',[ln_mu_start ln_sigma_start]);
                ft_ln = fittype('logncdf(x,mu,sigma,''upper'')','options',fo_ln);
                [cf_ln,~] = fit(X,ccdf,ft_ln);
                thispara_on(i,:) = coeffvalues(cf_ln);
            catch
                thispara_on(i,:) = [0 1];
            end
        elseif strcmp(ondist,'MittagLeffler')
            ml_beta_start = 0.5;
            ml_gamma_start = 0.5;
            try
                fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[ml_beta_start ml_gamma_start]);
                ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
                [cf_ml,~] = fit(X,ccdf,ft_ml);
                thispara_on(i,:) = coeffvalues(cf_ml);
            catch
                thispara_on(i,:) = [ml_beta_start ml_gamma_start];
            end
        elseif strcmp(ondist,'GenPareto')
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
            try
                fo_gp = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf -inf 0],'StartPoint',[gp_k_start gp_sigma_start gp_theta_start]);
                ft_gp = fittype('gpcdf(x,k,sigma,theta,''upper'')','options',fo_gp);
                [cf_gp,~] = fit(X,ccdf,ft_gp);
                thispara_on(i,:) = coeffvalues(cf_gp);
            catch
                thispara_on(i,:) = [gp_k_start gp_sigma_start gp_theta_start];
            end
        elseif strcmp(ondist,'Weibull')
            wb_a_start = 0.5;
            wb_b_start = 0.5;
            try
                fo_wb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'StartPoint',[wb_a_start wb_b_start]);
                ft_wb = fittype('wblcdf(x,a,b,''upper'')','options',fo_wb);
                [cf_wb,~] = fit(X,ccdf,ft_wb);
                thispara_on(i,:) = coeffvalues(cf_wb);
            catch
                thispara_on(i,:) = [wb_a_start wb_b_start];
            end
        end
        
    end
    
    %Fitting for off-periods
    fields_off = fieldnames(OffPeriods);
    
    thispara_off = zeros(numel(fields_off),paracount_off);
    
    
    for i = 1:numel(fields_off)
        thisoff = OnPeriods.(fields_off{i});
        [F,X] = ecdf(thisoff);
        ccdf = 1-F;
        M1 = mean(thisoff);
        M2 = mean(thisoff.^2);
        M3 = mean(thisoff.^3);
        
        if strcmp(offdist,'Exponential')
            ex_lambda_start = M1;
            if ex_lambda_start<=0 || isnan(ex_lambda_start) || isinf(ex_lambda_start)
                ex_lambda_start = 0.1;
            end
            try
                fo_ex = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[ex_lambda_start]);
                ft_ex = fittype('expcdf(x,lambda,''upper'')','options',fo_ex);
                [cf_ex,~] = fit(X,ccdf,ft_ex);
                thispara_off(i,:) = coeffvalues(cf_ex);
            catch
                thispara_off(i,:) = [ex_lambda_start];
            end
        elseif strcmp(offdist,'Gamma')
            gm_b_start = (M2/M1)-M1;
            gm_a_start = M1/gm_b_start;
            if gm_b_start<=0 || isnan(gm_b_start) || isinf(gm_b_start)
                gm_b_start = 0.1;
            end
            if gm_a_start<=0 || isnan(gm_a_start) || isinf(gm_a_start)
                gm_a_start = 0.1;
            end
            try
                fo_gm = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[inf inf],'StartPoint',[gm_a_start gm_b_start]);
                ft_gm = fittype('gamcdf(x,a,b,''upper'')','options',fo_gm);
                [cf_gm,~] = fit(X,ccdf,ft_gm);
                thispara_off(i,:) = coeffvalues(cf_gm);
            catch
                thispara_off(i,:) = [gm_a_start gm_b_start];
            end
        elseif strcmp(offdist,'Rayleigh')
            rl_sigma_start = M1*sqrt(2/pi);
            if rl_sigma_start<=0 || isnan(rl_sigma_start) || isinf(rl_sigma_start)
                rl_sigma_start = 0.1;
            end
            try
                fo_rl = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0],'Upper',[inf],'StartPoint',[rl_sigma_start]);
                ft_rl = fittype('raylcdf(x,sigma,''upper'')','options',fo_rl);
                [cf_rl,~] = fit(X,ccdf,ft_rl);
                thispara_off(i,:) = coeffvalues(cf_rl);
            catch
                thispara_off(i,:) = [rl_sigma_start];
            end
        elseif strcmp(offdist,'LogNormal')
            ln_sigma_start = sqrt(log(M2*(M1^-2)));
            if ln_sigma_start<=0 || isnan(ln_sigma_start) || isinf(ln_sigma_start)
                ln_sigma_start = 0.1;
            end
            ln_mu_start = log(M1)-(0.5*ln_sigma_start^2);
            if isnan(ln_mu_start) || isinf(ln_mu_start)
                ln_mu_start = 0;
            end
            try
                fo_ln = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf 0],'Upper',[inf inf],'StartPoint',[ln_mu_start ln_sigma_start]);
                ft_ln = fittype('logncdf(x,mu,sigma,''upper'')','options',fo_ln);
                [cf_ln,~] = fit(X,ccdf,ft_ln);
                thispara_off(i,:) = coeffvalues(cf_ln);
            catch
                thispara_off(i,:) = [0 1];
            end
        elseif strcmp(offdist,'MittagLeffler')
            ml_beta_start = 0.5;
            ml_gamma_start = 0.5;
            try
                fo_ml = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'Upper',[1 1],'StartPoint',[ml_beta_start ml_gamma_start]);
                ft_ml = fittype('mlf(beta,1,-gamma*x.^beta,6)','options',fo_ml);
                [cf_ml,~] = fit(X,ccdf,ft_ml);
                thispara_off(i,:) = coeffvalues(cf_ml);
            catch
                thispara_off(i,:) = [ml_beta_start ml_gamma_start];
            end
        elseif strcmp(offdist,'GenPareto')
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
            try
                fo_gp = fitoptions('Method', 'NonlinearLeastSquares','Lower',[-inf -inf 0],'StartPoint',[gp_k_start gp_sigma_start gp_theta_start]);
                ft_gp = fittype('gpcdf(x,k,sigma,theta,''upper'')','options',fo_gp);
                [cf_gp,~] = fit(X,ccdf,ft_gp);
                thispara_off(i,:) = coeffvalues(cf_gp);
            catch
                thispara_off(i,:) = [gp_k_start gp_sigma_start gp_theta_start];
            end
        elseif strcmp(offdist,'Weibull')
            wb_a_start = 0.5;
            wb_b_start = 0.5;
            try
                fo_wb = fitoptions('Method', 'NonlinearLeastSquares','Lower',[0 0],'StartPoint',[wb_a_start wb_b_start]);
                ft_wb = fittype('wblcdf(x,a,b,''upper'')','options',fo_wb);
                [cf_wb,~] = fit(X,ccdf,ft_wb);
                thispara_off(i,:) = coeffvalues(cf_wb);
            catch
                thispara_off(i,:) = [wb_a_start wb_b_start];
            end
        end
        
    end
    para_on = [para_on;thispara_on];
    para_off = [para_off;thispara_off];
end

paraon = struct();
for i=1:paracount_on
    thisdata = para_on(:,i)';
    [F,X] = ecdf(thisdata);
    ccdf = 1-F;
    M1 = mean(thisdata);
    M2 = mean(thisdata.^2);
    M3 = mean(thisdata.^3);
    
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
    
    sorteddata = sort(thisdata)';
    
    z_ex = expcdf(sorteddata,ex_lambda);
    z_gm = gamcdf(sorteddata,gm_a,gm_b);
    z_rl = raylcdf(sorteddata,rl_sigma);
    z_ln = logncdf(sorteddata,ln_mu,ln_sigma);
    z_ml = ones(length(sorteddata),1)-mlf(ml_beta,1,-ml_gamma*sorteddata.^ml_beta,6);
    z_gp = gpcdf(sorteddata,gp_k,gp_sigma,gp_theta);
    z_wb = wblcdf(sorteddata,wb_a,wb_b);
    
    zp_ex = exppdf(sorteddata,ex_lambda);
    zp_gm = gampdf(sorteddata,gm_a,gm_b);
    zp_rl = raylpdf(sorteddata,rl_sigma);
    zp_ln = lognpdf(sorteddata,ln_mu,ln_sigma);
    zp_ml = (-ml_beta./sorteddata).*mlf(ml_beta,1,-ml_gamma*sorteddata.^ml_beta,6);
    zp_gp = gppdf(sorteddata,gp_k,gp_sigma,gp_theta);
    zp_wb = wblpdf(sorteddata,wb_a,wb_b);
    
    stats_ex = testStatistics(sorteddata,z_ex,zp_ex,20);
    stats_gm = testStatistics(sorteddata,z_gm,zp_gm,20);
    stats_rl = testStatistics(sorteddata,z_rl,zp_rl,20);
    stats_ln = testStatistics(sorteddata,z_ln,zp_ln,20);
    stats_ml = testStatistics(sorteddata,z_ml,zp_ml,20);
    stats_gp = testStatistics(sorteddata,z_gp,zp_gp,20);
    stats_wb = testStatistics(sorteddata,z_wb,zp_wb,20);
    
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
    name = sprintf('Para%d', i);
    paraon.(name) = thisStructure;
    
    ccdf_ex = expcdf(X,ex_lambda,'upper');
    ccdf_gm = gamcdf(X,gm_a,gm_b,'upper');
    ccdf_rl = raylcdf(X,rl_sigma,'upper');
    ccdf_ln = logncdf(X,ln_mu,ln_sigma,'upper');
    ccdf_ml = mlf(ml_beta,1,-ml_gamma*X.^ml_beta,6);
    ccdf_gp = gpcdf(X,gp_k,gp_sigma,gp_theta,'upper');
    ccdf_wb = wblcdf(X,wb_a,wb_b,'upper');
    
    fig = figure();
    hold on
    plot(X,ccdf,'o');
    plot(X,ccdf_ex);
    plot(X,ccdf_gm);
    plot(X,ccdf_rl);
    plot(X,ccdf_ln);
    plot(X,ccdf_ml);
    plot(X,ccdf_gp);
    plot(X,ccdf_wb);
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    xlabel('Parameter Values');
    ylabel('CCDF');
    axis([-inf,inf,1E-10,1E0]);
    legend('Data','Exponential','Gamma','Rayleigh','Log-Normal','Mittag-Leffler','Gen. Pareto','Weibull','Location','southwest');
    imagefilename = [dir_ref,'/on_para',num2str(i),'.png'];
    print(imagefilename,'-dpng')
    close(fig);
end

paraoff = struct();
for i=1:paracount_off
    thisdata = para_off(:,i)';
    [F,X] = ecdf(thisdata);
    ccdf = 1-F;
    M1 = mean(thisdata);
    M2 = mean(thisdata.^2);
    M3 = mean(thisdata.^3);
    
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
    
    sorteddata = sort(thisdata)';
    
    z_ex = expcdf(sorteddata,ex_lambda);
    z_gm = gamcdf(sorteddata,gm_a,gm_b);
    z_rl = raylcdf(sorteddata,rl_sigma);
    z_ln = logncdf(sorteddata,ln_mu,ln_sigma);
    z_ml = ones(length(sorteddata),1)-mlf(ml_beta,1,-ml_gamma*sorteddata.^ml_beta,6);
    z_gp = gpcdf(sorteddata,gp_k,gp_sigma,gp_theta);
    z_wb = wblcdf(sorteddata,wb_a,wb_b);
    
    zp_ex = exppdf(sorteddata,ex_lambda);
    zp_gm = gampdf(sorteddata,gm_a,gm_b);
    zp_rl = raylpdf(sorteddata,rl_sigma);
    zp_ln = lognpdf(sorteddata,ln_mu,ln_sigma);
    zp_ml = (-ml_beta./sorteddata).*mlf(ml_beta,1,-ml_gamma*sorteddata.^ml_beta,6);
    zp_gp = gppdf(sorteddata,gp_k,gp_sigma,gp_theta);
    zp_wb = wblpdf(sorteddata,wb_a,wb_b);
    
    stats_ex = testStatistics(sorteddata,z_ex,zp_ex,20);
    stats_gm = testStatistics(sorteddata,z_gm,zp_gm,20);
    stats_rl = testStatistics(sorteddata,z_rl,zp_rl,20);
    stats_ln = testStatistics(sorteddata,z_ln,zp_ln,20);
    stats_ml = testStatistics(sorteddata,z_ml,zp_ml,20);
    stats_gp = testStatistics(sorteddata,z_gp,zp_gp,20);
    stats_wb = testStatistics(sorteddata,z_wb,zp_wb,20);
    
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
    name = sprintf('Para%d', i);
    paraoff.(name) = thisStructure;
        
    ccdf_ex = expcdf(X,ex_lambda,'upper');
    ccdf_gm = gamcdf(X,gm_a,gm_b,'upper');
    ccdf_rl = raylcdf(X,rl_sigma,'upper');
    ccdf_ln = logncdf(X,ln_mu,ln_sigma,'upper');
    ccdf_ml = mlf(ml_beta,1,-ml_gamma*X.^ml_beta,6);
    ccdf_gp = gpcdf(X,gp_k,gp_sigma,gp_theta,'upper');
    ccdf_wb = wblcdf(X,wb_a,wb_b,'upper');
    
    fig = figure();
    hold on
    plot(X,ccdf,'o');
    plot(X,ccdf_ex);
    plot(X,ccdf_gm);
    plot(X,ccdf_rl);
    plot(X,ccdf_ln);
    plot(X,ccdf_ml);
    plot(X,ccdf_gp);
    plot(X,ccdf_wb);
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    xlabel('Parameter Values');
    ylabel('CCDF');
    axis([-inf,inf,1E-10,1E0]);
    legend('Data','Exponential','Gamma','Rayleigh','Log-Normal','Mittag-Leffler','Gen. Pareto','Weibull','Location','southwest');
    imagefilename = [dir_ref,'/off_para',num2str(i),'.png'];
    print(imagefilename,'-dpng')
    close(fig);
end

parameterAnalysisLaTeX(paraon,dir_ref,'OnDistributions_paras',ondist);
parameterAnalysisLaTeX(paraoff,dir_ref,'OffDistributions_paras',offdist);

end