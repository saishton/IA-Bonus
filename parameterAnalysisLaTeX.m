function [] = parameterAnalysisLaTeX(structure,dir_ref,name,dist)

paraNames = struct();
if strcmp(dist,'Exponential')
    paracount = 1;
    paraNames.p1 = 'Scale';
elseif strcmp(dist,'Gamma')
    paracount = 2;
    paraNames.p1 = 'Shape';
    paraNames.p2 = 'Scale';
elseif strcmp(dist,'Rayleigh')
    paracount = 1;
    paraNames.p1 = 'Scale';
elseif strcmp(dist,'LogNormal')
    paracount = 2;
    paraNames.p1 = 'Location';
    paraNames.p2 = 'Scale';
elseif strcmp(dist,'MittagLeffler')
    paracount = 2;
    paraNames.p1 = 'Stability';
    paraNames.p2 = 'Scale';
elseif strcmp(dist,'GenPareto')
    paracount = 3;
    paraNames.p1 = 'Shape';
    paraNames.p2 = 'Scale';
    paraNames.p3 = 'Location';
elseif strcmp(dist,'Weibull')
    paracount = 2;
    paraNames.p1 = 'Scale';
    paraNames.p2 = 'Shape';
end

cleanName = strrep(name, '.', '');
cleanName = strrep(cleanName, ' ', '_');

parfor i=1:paracount
    filename = ['LaTexTables_Distributions_',cleanName,'_parameter_',num2str(i),'.txt'];
    filepath = [dir_ref,'/',filename];
    
    pararef = sprintf('p%d', i);
    paraname = paraNames.(pararef);
    strname = sprintf('Para%d', i);
    paraStr = structure.(strname);
    EX = paraStr.Exponential;
    GM = paraStr.Gamma;
    RL = paraStr.Rayleigh;
    LN = paraStr.LogNormal;
    ML = paraStr.MittagLeffler;
    GP = paraStr.GenPareto;
    WB = paraStr.Weibull;
        
    lines = 18;
    tobuild = cell(lines,1);
    
    tobuild{01} = ['\subsection{',dist,': ',paraname,'}'];
    tobuild{02} = '\begin{tabular}{|c||c|c||c|c|c|c|c|c|c|} \hline';
    tobuild{03} = '\multirow{2}{*}{Distribution} & \multicolumn{2}{|c||}{\multirow{2}{*}{Parameters}} & \multicolumn{7}{|c|}{Statistics} \\ \cline{4-10}';
    tobuild{04} = ' & \multicolumn{2}{|c||}{} & KolD & CvM & Kuiper & Watson & And-Dar & Kull-Lei & Jen-Sha \\ \hline';
    tobuild{05} = ['Exponential & Scale & ',num2matlabstr(EX.Parameters.Scale),' & ', latexstats(EX.Statistics,1),' \\ \hline'];
    tobuild{06} = ['\multirow{2}{*}{Gamma} & Shape & ',num2matlabstr(GM.Parameters.Shape),' & ', latexstats(GM.Statistics,2),' \\ \cline{2-3}'];
    tobuild{07} = [' & Scale & ',num2matlabstr(GM.Parameters.Scale),' & & & & & & & \\ \hline'];
    tobuild{08} = ['Rayleigh & Scale & ',num2matlabstr(RL.Parameters.Scale),' & ', latexstats(RL.Statistics,1),' \\ \hline'];
    tobuild{09} = ['\multirow{2}{*}{Log-Normal} & Location & ',num2matlabstr(LN.Parameters.Location),' & ', latexstats(LN.Statistics,2),' \\ \cline{2-3}'];
    tobuild{10} = [' & Scale & ',num2matlabstr(LN.Parameters.Scale),' & & & & & & & \\ \hline'];
    tobuild{11} = ['\multirow{2}{*}{Mittag-Leffler} & Stability & ',num2matlabstr(ML.Parameters.Stability),' & ', latexstats(ML.Statistics,2),' \\ \cline{2-3}'];
    tobuild{12} = [' & Scale & ',num2matlabstr(ML.Parameters.Scale),' & & & & & & & \\ \hline'];
    tobuild{13} = ['\multirow{3}{*}{Gen. Pareto} & Shape & ',num2matlabstr(GP.Parameters.Shape),' & ', latexstats(GP.Statistics,3),' \\ \cline{2-3}'];
    tobuild{14} = [' & Scale & ',num2matlabstr(GP.Parameters.Scale),' & & & & & & & \\ \cline{2-3}'];
    tobuild{15} = [' & Location & ',num2matlabstr(GP.Parameters.Location),' & & & & & & & \\ \hline'];
    tobuild{16} = ['\multirow{2}{*}{Weibull} & Scale & ',num2matlabstr(WB.Parameters.Scale),' & ', latexstats(WB.Statistics,2),' \\ \cline{2-3}'];
    tobuild{17} = [' & Shape & ',num2matlabstr(WB.Parameters.Shape),' & & & & & & & \\ \hline'];
    tobuild{18} = '\end{tabular}';
      
    
    fileID = fopen(filepath,'w');
    fprintf(fileID,'%s\r\n',tobuild{:});
    fclose(fileID);
end

end