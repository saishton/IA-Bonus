function [statstr] = latexstats(stats,multirow)

if multirow == 1
    statstr = [num2matlabstr(stats.Kolmogorov_D),' & ',num2matlabstr(stats.Cramer_von_Mises),' & ',num2matlabstr(stats.Kuiper),' & ',num2matlabstr(stats.Watson),' & ',num2matlabstr(stats.Anderson_Darling),' & ',num2matlabstr(stats.Kullback_Leibler),' & ',num2matlabstr(stats.Jensen_Shannon)];
else
    statstr = ['\multirow{',num2str(multirow),'}{*}{',num2matlabstr(stats.Kolmogorov_D),'} & \multirow{',num2str(multirow),'}{*}{',num2matlabstr(stats.Cramer_von_Mises),'} & \multirow{',num2str(multirow),'}{*}{',num2matlabstr(stats.Kuiper),'} & \multirow{',num2str(multirow),'}{*}{',num2matlabstr(stats.Watson),'} & \multirow{',num2str(multirow),'}{*}{',num2matlabstr(stats.Anderson_Darling),'} & \multirow{',num2str(multirow),'}{*}{',num2matlabstr(stats.Kullback_Leibler),'} & \multirow{',num2str(multirow),'}{*}{',num2matlabstr(stats.Jensen_Shannon),'}'];
end
end