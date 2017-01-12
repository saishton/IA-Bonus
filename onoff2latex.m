function [] = onoff2latex(analysis,dir_ref,name)

[ssum,avrank] = statSum(analysis);

cleanName = strrep(name, '.', '');
cleanName = strrep(cleanName, ' ', '_');

filename = ['LaTexTables_Distributions_',cleanName,'_stats.txt'];
filepath = [dir_ref,'/',filename];

filename2 = ['LaTexTables_Distributions_',cleanName,'_comps.txt'];
filepath2 = [dir_ref,'/',filename2];

names = fieldnames(analysis);
samples = numel(names);
prebonus = 3;
postbonus = 3;
lines = samples+prebonus+postbonus;

stat_string = '& \rotatebox[origin=c]{90}{Kolmogorov D} & \rotatebox[origin=c]{90}{Cramer Von Mises} & \rotatebox[origin=c]{90}{Kuiper} & \rotatebox[origin=c]{90}{Watson} & \rotatebox[origin=c]{90}{Anderson-Darling} & \rotatebox[origin=c]{90}{Kullback-Leibler} & \rotatebox[origin=c]{90}{Jensen-Shannon} \\ \hline';

EX_build = cell(lines,1);
GM_build = cell(lines,1);
RL_build = cell(lines,1);
LN_build = cell(lines,1);
ML_build = cell(lines,1);
GP_build = cell(lines,1);
WB_build = cell(lines,1);

EX_build{1} = '\subsection{Results: Exponential}';
GM_build{1} = '\subsection{Results: Gamma}';
RL_build{1} = '\subsection{Results: Rayleigh}';
LN_build{1} = '\subsection{Results: Log-Normal}';
ML_build{1} = '\subsection{Results: Mittag-Leffler}';
GP_build{1} = '\subsection{Results: Generalised Pareto}';
WB_build{1} = '\subsection{Results: Weibull}';

EX_build{2} = '\begin{tabular}{|c||c||c|c|c|c|c|c|c|} \hline';
GM_build{2} = '\begin{tabular}{|c||c|c||c|c|c|c|c|c|c|} \hline';
RL_build{2} = '\begin{tabular}{|c||c||c|c|c|c|c|c|c|} \hline';
LN_build{2} = '\begin{tabular}{|c||c|c||c|c|c|c|c|c|c|} \hline';
ML_build{2} = '\begin{tabular}{|c||c|c||c|c|c|c|c|c|c|} \hline';
GP_build{2} = '\begin{tabular}{|c||c|c|c||c|c|c|c|c|c|c|} \hline';
WB_build{2} = '\begin{tabular}{|c||c|c||c|c|c|c|c|c|c|} \hline';

EX_build{3} = ['\rotatebox[origin=c]{90}{Connection} & \rotatebox[origin=c]{90}{Scale} ',stat_string]; 
GM_build{3} = ['\rotatebox[origin=c]{90}{Connection} & \rotatebox[origin=c]{90}{Shape} & \rotatebox[origin=c]{90}{Scale} ',stat_string]; 
RL_build{3} = ['\rotatebox[origin=c]{90}{Connection} & \rotatebox[origin=c]{90}{Scale} ',stat_string]; 
LN_build{3} = ['\rotatebox[origin=c]{90}{Connection} & \rotatebox[origin=c]{90}{Location} & \rotatebox[origin=c]{90}{Scale} ',stat_string]; 
ML_build{3} = ['\rotatebox[origin=c]{90}{Connection} & \rotatebox[origin=c]{90}{Stability} & \rotatebox[origin=c]{90}{Scale} ',stat_string]; 
GP_build{3} = ['\rotatebox[origin=c]{90}{Connection} & \rotatebox[origin=c]{90}{Shape} & \rotatebox[origin=c]{90}{Scale} & \rotatebox[origin=c]{90}{Location} ',stat_string]; 
WB_build{3} = ['\rotatebox[origin=c]{90}{Connection} & \rotatebox[origin=c]{90}{Scale} & \rotatebox[origin=c]{90}{Shape} ',stat_string]; 

parfor i=1:samples
    
    con = strrep(names{i}, '_', '\leftrightarrow');
    con = strrep(con, 'n', '');
    con = ['$',con,'$'];
    
    thisdata = analysis.(names{i});
    
    thisEX = thisdata.Exponential;
    thisGM = thisdata.Gamma;
    thisRL = thisdata.Rayleigh;
    thisLN = thisdata.LogNormal;
    thisML = thisdata.MittagLeffler;
    thisGP = thisdata.GenPareto;
    thisWB = thisdata.Weibull;
    
    thisEXpara = thisEX.Parameters;
    thisGMpara = thisGM.Parameters;
    thisRLpara = thisRL.Parameters;
    thisLNpara = thisLN.Parameters;
    thisMLpara = thisML.Parameters;
    thisGPpara = thisGP.Parameters;
    thisWBpara = thisWB.Parameters;
    
    thisEXstat = thisEX.Statistics;
    thisGMstat = thisGM.Statistics;
    thisRLstat = thisRL.Statistics;
    thisLNstat = thisLN.Statistics;
    thisMLstat = thisML.Statistics;
    thisGPstat = thisGP.Statistics;
    thisWBstat = thisWB.Statistics;
    
    EXstatStr = [' & ',num2matlabstr(thisEXstat.Kolmogorov_D),' & ',num2matlabstr(thisEXstat.Cramer_von_Mises),' & ',num2matlabstr(thisEXstat.Kuiper),' & ',num2matlabstr(thisEXstat.Watson),' & ',num2matlabstr(thisEXstat.Anderson_Darling),' & ',num2matlabstr(thisEXstat.Kullback_Leibler),' & ',num2matlabstr(thisEXstat.Jensen_Shannon)];
    GMstatStr = [' & ',num2matlabstr(thisGMstat.Kolmogorov_D),' & ',num2matlabstr(thisGMstat.Cramer_von_Mises),' & ',num2matlabstr(thisGMstat.Kuiper),' & ',num2matlabstr(thisGMstat.Watson),' & ',num2matlabstr(thisGMstat.Anderson_Darling),' & ',num2matlabstr(thisGMstat.Kullback_Leibler),' & ',num2matlabstr(thisGMstat.Jensen_Shannon)];
    RLstatStr = [' & ',num2matlabstr(thisRLstat.Kolmogorov_D),' & ',num2matlabstr(thisRLstat.Cramer_von_Mises),' & ',num2matlabstr(thisRLstat.Kuiper),' & ',num2matlabstr(thisRLstat.Watson),' & ',num2matlabstr(thisRLstat.Anderson_Darling),' & ',num2matlabstr(thisRLstat.Kullback_Leibler),' & ',num2matlabstr(thisRLstat.Jensen_Shannon)];
    LNstatStr = [' & ',num2matlabstr(thisLNstat.Kolmogorov_D),' & ',num2matlabstr(thisLNstat.Cramer_von_Mises),' & ',num2matlabstr(thisLNstat.Kuiper),' & ',num2matlabstr(thisLNstat.Watson),' & ',num2matlabstr(thisLNstat.Anderson_Darling),' & ',num2matlabstr(thisLNstat.Kullback_Leibler),' & ',num2matlabstr(thisLNstat.Jensen_Shannon)];
    MLstatStr = [' & ',num2matlabstr(thisMLstat.Kolmogorov_D),' & ',num2matlabstr(thisMLstat.Cramer_von_Mises),' & ',num2matlabstr(thisMLstat.Kuiper),' & ',num2matlabstr(thisMLstat.Watson),' & ',num2matlabstr(thisMLstat.Anderson_Darling),' & ',num2matlabstr(thisMLstat.Kullback_Leibler),' & ',num2matlabstr(thisMLstat.Jensen_Shannon)];
    GPstatStr = [' & ',num2matlabstr(thisGPstat.Kolmogorov_D),' & ',num2matlabstr(thisGPstat.Cramer_von_Mises),' & ',num2matlabstr(thisGPstat.Kuiper),' & ',num2matlabstr(thisGPstat.Watson),' & ',num2matlabstr(thisGPstat.Anderson_Darling),' & ',num2matlabstr(thisGPstat.Kullback_Leibler),' & ',num2matlabstr(thisGPstat.Jensen_Shannon)];
    WBstatStr = [' & ',num2matlabstr(thisWBstat.Kolmogorov_D),' & ',num2matlabstr(thisWBstat.Cramer_von_Mises),' & ',num2matlabstr(thisWBstat.Kuiper),' & ',num2matlabstr(thisWBstat.Watson),' & ',num2matlabstr(thisWBstat.Anderson_Darling),' & ',num2matlabstr(thisWBstat.Kullback_Leibler),' & ',num2matlabstr(thisWBstat.Jensen_Shannon)];
    
    EX_build{prebonus+i} = [con,' & ',num2matlabstr(thisEXpara.Scale),EXstatStr,'\\ \hline'];
    GM_build{prebonus+i} = [con,' & ',num2matlabstr(thisGMpara.Shape),' & ',num2matlabstr(thisGMpara.Scale),GMstatStr,'\\ \hline'];
    RL_build{prebonus+i} = [con,' & ',num2matlabstr(thisRLpara.Scale),RLstatStr,'\\ \hline'];
    LN_build{prebonus+i} = [con,' & ',num2matlabstr(thisLNpara.Location),' & ',num2matlabstr(thisLNpara.Scale),LNstatStr,'\\ \hline'];
    ML_build{prebonus+i} = [con,' & ',num2matlabstr(thisMLpara.Stability),' & ',num2matlabstr(thisMLpara.Scale),MLstatStr,'\\ \hline'];
    GP_build{prebonus+i} = [con,' & ',num2matlabstr(thisGPpara.Shape),' & ',num2matlabstr(thisGPpara.Scale),' & ',num2matlabstr(thisGPpara.Location),GPstatStr,'\\ \hline'];
    WB_build{prebonus+i} = [con,' & ',num2matlabstr(thisWBpara.Scale),' & ',num2matlabstr(thisWBpara.Shape),WBstatStr,'\\ \hline'];
    
end

EX_build{samples+prebonus+1} = ['\textbf{Sum} & & ',num2matlabstr(ssum(1,1)),' & ',num2matlabstr(ssum(1,2)),' & ',num2matlabstr(ssum(1,3)),' & ',num2matlabstr(ssum(1,4)),' & ',num2matlabstr(ssum(1,5)),' & ',num2matlabstr(ssum(1,6)),' & ',num2matlabstr(ssum(1,7)),'\\ \hline'];
GM_build{samples+prebonus+1} = ['\textbf{Sum} & & & ',num2matlabstr(ssum(2,1)),' & ',num2matlabstr(ssum(2,2)),' & ',num2matlabstr(ssum(2,3)),' & ',num2matlabstr(ssum(2,4)),' & ',num2matlabstr(ssum(2,5)),' & ',num2matlabstr(ssum(2,6)),' & ',num2matlabstr(ssum(2,7)),'\\ \hline'];
RL_build{samples+prebonus+1} = ['\textbf{Sum} & & ',num2matlabstr(ssum(3,1)),' & ',num2matlabstr(ssum(3,2)),' & ',num2matlabstr(ssum(3,3)),' & ',num2matlabstr(ssum(3,4)),' & ',num2matlabstr(ssum(3,5)),' & ',num2matlabstr(ssum(3,6)),' & ',num2matlabstr(ssum(3,7)),'\\ \hline'];
LN_build{samples+prebonus+1} = ['\textbf{Sum} & & & ',num2matlabstr(ssum(4,1)),' & ',num2matlabstr(ssum(4,2)),' & ',num2matlabstr(ssum(4,3)),' & ',num2matlabstr(ssum(4,4)),' & ',num2matlabstr(ssum(4,5)),' & ',num2matlabstr(ssum(4,6)),' & ',num2matlabstr(ssum(4,7)),'\\ \hline'];
ML_build{samples+prebonus+1} = ['\textbf{Sum} & & & ',num2matlabstr(ssum(5,1)),' & ',num2matlabstr(ssum(5,2)),' & ',num2matlabstr(ssum(5,3)),' & ',num2matlabstr(ssum(5,4)),' & ',num2matlabstr(ssum(5,5)),' & ',num2matlabstr(ssum(5,6)),' & ',num2matlabstr(ssum(5,7)),'\\ \hline'];
GP_build{samples+prebonus+1} = ['\textbf{Sum} & & & & ',num2matlabstr(ssum(6,1)),' & ',num2matlabstr(ssum(6,2)),' & ',num2matlabstr(ssum(6,3)),' & ',num2matlabstr(ssum(6,4)),' & ',num2matlabstr(ssum(6,5)),' & ',num2matlabstr(ssum(6,6)),' & ',num2matlabstr(ssum(6,7)),'\\ \hline'];
WB_build{samples+prebonus+1} = ['\textbf{Sum} & & & ',num2matlabstr(ssum(7,1)),' & ',num2matlabstr(ssum(7,2)),' & ',num2matlabstr(ssum(7,3)),' & ',num2matlabstr(ssum(7,4)),' & ',num2matlabstr(ssum(7,5)),' & ',num2matlabstr(ssum(7,6)),' & ',num2matlabstr(ssum(7,7)),'\\ \hline'];

EX_build{samples+prebonus+2} = '\end{tabular}';
GM_build{samples+prebonus+2} = '\end{tabular}';
RL_build{samples+prebonus+2} = '\end{tabular}';
LN_build{samples+prebonus+2} = '\end{tabular}';
ML_build{samples+prebonus+2} = '\end{tabular}';
GP_build{samples+prebonus+2} = '\end{tabular}';
WB_build{samples+prebonus+2} = '\end{tabular}';

EX_build{samples+prebonus+3} = '\newpage';
GM_build{samples+prebonus+3} = '\newpage';
RL_build{samples+prebonus+3} = '\newpage';
LN_build{samples+prebonus+3} = '\newpage';
ML_build{samples+prebonus+3} = '\newpage';
GP_build{samples+prebonus+3} = '\newpage';
WB_build{samples+prebonus+3} = '\newpage';

comp_build = cell(22,1);
comp_build{1} = '\subsection{Results: Comparison of Stats}';
comp_build{2} = '\begin{tabular}{|c||c|c|c|c|c|c|c|} \hline';
comp_build{3} = ['\rotatebox[origin=c]{90}{Distribution} ',stat_string];
comp_build{4} = ['\textbf{Exponential} & ',num2matlabstr(ssum(1,1)),' & ',num2matlabstr(ssum(1,2)),' & ',num2matlabstr(ssum(1,3)),' & ',num2matlabstr(ssum(1,4)),' & ',num2matlabstr(ssum(1,5)),' & ',num2matlabstr(ssum(1,6)),' & ',num2matlabstr(ssum(1,7)),'\\ \hline'];
comp_build{5} = ['\textbf{Gamma} & ',num2matlabstr(ssum(2,1)),' & ',num2matlabstr(ssum(2,2)),' & ',num2matlabstr(ssum(2,3)),' & ',num2matlabstr(ssum(2,4)),' & ',num2matlabstr(ssum(2,5)),' & ',num2matlabstr(ssum(2,6)),' & ',num2matlabstr(ssum(2,7)),'\\ \hline'];
comp_build{6} = ['\textbf{Rayleigh} & ',num2matlabstr(ssum(3,1)),' & ',num2matlabstr(ssum(3,2)),' & ',num2matlabstr(ssum(3,3)),' & ',num2matlabstr(ssum(3,4)),' & ',num2matlabstr(ssum(3,5)),' & ',num2matlabstr(ssum(3,6)),' & ',num2matlabstr(ssum(3,7)),'\\ \hline'];
comp_build{7} = ['\textbf{Log-Normal} & ',num2matlabstr(ssum(4,1)),' & ',num2matlabstr(ssum(4,2)),' & ',num2matlabstr(ssum(4,3)),' & ',num2matlabstr(ssum(4,4)),' & ',num2matlabstr(ssum(4,5)),' & ',num2matlabstr(ssum(4,6)),' & ',num2matlabstr(ssum(4,7)),'\\ \hline'];
comp_build{8} = ['\textbf{Mittag-Leffler} & ',num2matlabstr(ssum(5,1)),' & ',num2matlabstr(ssum(5,2)),' & ',num2matlabstr(ssum(5,3)),' & ',num2matlabstr(ssum(5,4)),' & ',num2matlabstr(ssum(5,5)),' & ',num2matlabstr(ssum(5,6)),' & ',num2matlabstr(ssum(5,7)),'\\ \hline'];
comp_build{9} = ['\textbf{Generalised Pareto} & ',num2matlabstr(ssum(6,1)),' & ',num2matlabstr(ssum(6,2)),' & ',num2matlabstr(ssum(6,3)),' & ',num2matlabstr(ssum(6,4)),' & ',num2matlabstr(ssum(6,5)),' & ',num2matlabstr(ssum(6,6)),' & ',num2matlabstr(ssum(6,7)),'\\ \hline'];
comp_build{10} = ['\textbf{Weibull} & ',num2matlabstr(ssum(7,1)),' & ',num2matlabstr(ssum(7,2)),' & ',num2matlabstr(ssum(7,3)),' & ',num2matlabstr(ssum(7,4)),' & ',num2matlabstr(ssum(7,5)),' & ',num2matlabstr(ssum(7,6)),' & ',num2matlabstr(ssum(7,7)),'\\ \hline'];
comp_build{11} = '\end{tabular}';
comp_build{12} = '\subsection{Results: Comparison of Average Rank}';
comp_build{13} = '\begin{tabular}{|c||c|c|c|c|c|c|c|} \hline';
comp_build{14} = ['\rotatebox[origin=c]{90}{Distribution} ',stat_string];
comp_build{15} = ['\textbf{Exponential} & ',num2matlabstr(avrank(1,1)),' & ',num2matlabstr(avrank(1,2)),' & ',num2matlabstr(avrank(1,3)),' & ',num2matlabstr(avrank(1,4)),' & ',num2matlabstr(avrank(1,5)),' & ',num2matlabstr(avrank(1,6)),' & ',num2matlabstr(avrank(1,7)),'\\ \hline'];
comp_build{16} = ['\textbf{Gamma} & ',num2matlabstr(avrank(2,1)),' & ',num2matlabstr(avrank(2,2)),' & ',num2matlabstr(avrank(2,3)),' & ',num2matlabstr(avrank(2,4)),' & ',num2matlabstr(avrank(2,5)),' & ',num2matlabstr(avrank(2,6)),' & ',num2matlabstr(avrank(2,7)),'\\ \hline'];
comp_build{17} = ['\textbf{Rayleigh} & ',num2matlabstr(avrank(3,1)),' & ',num2matlabstr(avrank(3,2)),' & ',num2matlabstr(avrank(3,3)),' & ',num2matlabstr(avrank(3,4)),' & ',num2matlabstr(avrank(3,5)),' & ',num2matlabstr(avrank(3,6)),' & ',num2matlabstr(avrank(3,7)),'\\ \hline'];
comp_build{18} = ['\textbf{Log-Normal} & ',num2matlabstr(avrank(4,1)),' & ',num2matlabstr(avrank(4,2)),' & ',num2matlabstr(avrank(4,3)),' & ',num2matlabstr(avrank(4,4)),' & ',num2matlabstr(avrank(4,5)),' & ',num2matlabstr(avrank(4,6)),' & ',num2matlabstr(avrank(4,7)),'\\ \hline'];
comp_build{19} = ['\textbf{Mittag-Leffler} & ',num2matlabstr(avrank(5,1)),' & ',num2matlabstr(avrank(5,2)),' & ',num2matlabstr(avrank(5,3)),' & ',num2matlabstr(avrank(5,4)),' & ',num2matlabstr(avrank(5,5)),' & ',num2matlabstr(avrank(5,6)),' & ',num2matlabstr(avrank(5,7)),'\\ \hline'];
comp_build{20} = ['\textbf{Generalised Pareto} & ',num2matlabstr(avrank(6,1)),' & ',num2matlabstr(avrank(6,2)),' & ',num2matlabstr(avrank(6,3)),' & ',num2matlabstr(avrank(6,4)),' & ',num2matlabstr(avrank(6,5)),' & ',num2matlabstr(avrank(6,6)),' & ',num2matlabstr(avrank(6,7)),'\\ \hline'];
comp_build{21} = ['\textbf{Weibull} & ',num2matlabstr(avrank(7,1)),' & ',num2matlabstr(avrank(7,2)),' & ',num2matlabstr(avrank(7,3)),' & ',num2matlabstr(avrank(7,4)),' & ',num2matlabstr(avrank(7,5)),' & ',num2matlabstr(avrank(7,6)),' & ',num2matlabstr(avrank(7,7)),'\\ \hline'];
comp_build{22} = '\end{tabular}';

tobuild = [EX_build; GM_build; RL_build; LN_build; ML_build; GP_build; WB_build];

fileID = fopen(filepath,'w');
fprintf(fileID,'%s\r\n',tobuild{:});
fclose(fileID);

fileID2 = fopen(filepath2,'w');
fprintf(fileID2,'%s\r\n',comp_build{:});
fclose(fileID2);
end

