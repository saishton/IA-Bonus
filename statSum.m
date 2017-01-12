function [ssum,avRank] = statSum(analysis)

ssum = zeros(7,7);
rsum = zeros(7,7);

names = fieldnames(analysis);
samples = numel(names);

for i=1:samples
    thisdata = analysis.(names{i});
    
    thisEX = thisdata.Exponential.Statistics;
    thisGM = thisdata.Gamma.Statistics;
    thisRL = thisdata.Rayleigh.Statistics;
    thisLN = thisdata.LogNormal.Statistics;
    thisML = thisdata.MittagLeffler.Statistics;
    thisGP = thisdata.GenPareto.Statistics;
    thisWB = thisdata.Weibull.Statistics;
    
    EXstats = [thisEX.Kolmogorov_D,thisEX.Cramer_von_Mises,thisEX.Kuiper,thisEX.Watson,thisEX.Anderson_Darling,thisEX.Kullback_Leibler,thisEX.Jensen_Shannon];
    GMstats = [thisGM.Kolmogorov_D,thisGM.Cramer_von_Mises,thisGM.Kuiper,thisGM.Watson,thisGM.Anderson_Darling,thisGM.Kullback_Leibler,thisGM.Jensen_Shannon];
    RLstats = [thisRL.Kolmogorov_D,thisRL.Cramer_von_Mises,thisRL.Kuiper,thisRL.Watson,thisRL.Anderson_Darling,thisRL.Kullback_Leibler,thisRL.Jensen_Shannon];
    LNstats = [thisLN.Kolmogorov_D,thisLN.Cramer_von_Mises,thisLN.Kuiper,thisLN.Watson,thisLN.Anderson_Darling,thisLN.Kullback_Leibler,thisLN.Jensen_Shannon];
    MLstats = [thisML.Kolmogorov_D,thisML.Cramer_von_Mises,thisML.Kuiper,thisML.Watson,thisML.Anderson_Darling,thisML.Kullback_Leibler,thisML.Jensen_Shannon];
    GPstats = [thisGP.Kolmogorov_D,thisGP.Cramer_von_Mises,thisGP.Kuiper,thisGP.Watson,thisGP.Anderson_Darling,thisGP.Kullback_Leibler,thisGP.Jensen_Shannon];
    WBstats = [thisWB.Kolmogorov_D,thisWB.Cramer_von_Mises,thisWB.Kuiper,thisWB.Watson,thisWB.Anderson_Darling,thisWB.Kullback_Leibler,thisWB.Jensen_Shannon];
    
    thisStats = [EXstats;GMstats;RLstats;LNstats;MLstats;GPstats;WBstats];
    ssum = ssum+thisStats;
    
    [~,~,rank_1] = unique(thisStats(:,1));
    [~,~,rank_2] = unique(thisStats(:,2));
    [~,~,rank_3] = unique(thisStats(:,3));
    [~,~,rank_4] = unique(thisStats(:,4));
    [~,~,rank_5] = unique(thisStats(:,5));
    [~,~,rank_6] = unique(thisStats(:,6));
    [~,~,rank_7] = unique(thisStats(:,7));
    
    rankMat = [rank_1,rank_2,rank_3,rank_4,rank_5,rank_6,rank_7];
    rsum = rsum+rankMat;
end

avRank = rsum./samples;

end
    