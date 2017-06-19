function [thisStructure] = fitdata(data)
[F,X] = ecdf(data);
ccdf = 1-F;
M1 = mean(data);
M2 = mean(data.^2);
M3 = mean(data.^3);

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

sortedon = sort(data)';

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

ccdf_ex = expcdf(X,ex_lambda,'upper');
ccdf_gm = gamcdf(X,gm_a,gm_b,'upper');
ccdf_rl = raylcdf(X,rl_sigma,'upper');
ccdf_ln = logncdf(X,ln_mu,ln_sigma,'upper');
ccdf_ml = mlf(ml_beta,1,-ml_gamma*X.^ml_beta,6);
ccdf_gp = gpcdf(X,gp_k,gp_sigma,gp_theta,'upper');
ccdf_wb = wblcdf(X,wb_a,wb_b,'upper');

graph = figure();
hold on
plot(X,ccdf,'x');
%plot(X,ccdf_ex);
%plot(X,ccdf_gm);
%plot(X,ccdf_rl);
plot(X,ccdf_ln);
%plot(X,ccdf_gp);
%plot(X,ccdf_wb);
set(gca,'XScale','log');
set(gca,'YScale','log');
ylim([min(ccdf(1:end-1)),max(ccdf)]);