function [AdjMat] = modelAM(initMat,runtime)

nodes = size(initMat,1);

ex1_mu = 3.2434;
EXpara1 = lognrnd(ex1_mu,sigma_for_mu_and_mean(30.552,ex1_mu),nodes);
LNpara1 = 6.3512*ones(nodes);
LNpara2 = 1.3688*ones(nodes);

AdjMat = zeros(nodes);

for i=1:nodes-1
    for j=i+1:nodes
        init = initMat(i,j);
        currenttime = 0;
        if init == 0
            thisoff = [0];
            thison = [];
            while currenttime<runtime
                thisoffduration = lognrnd(LNpara1(i,j),LNpara2(i,j));
                switch_on = currenttime+thisoffduration;
                thisonduration = exprnd(EXpara1(i,j));
                switch_off = switch_on+thisonduration;
                if switch_on<runtime
                    thison = [thison,switch_on];
                    if switch_off<runtime
                        thisoff = [thisoff,switch_off];
                    else
                        thisoff = [thisoff,runtime];
                    end
                else
                    thison = [thison,runtime];
                end
                currenttime = switch_off;
            end
        elseif init == 1
            thisoff = [];
            thison = [0];
            while currenttime<runtime
                thisonduration = exprnd(EXpara1(i,j));
                switch_off = currenttime+thisonduration;
                thisoffduration = lognrnd(LNpara1(i,j),LNpara2(i,j));
                switch_on = switch_off+thisoffduration;
                if switch_off<runtime
                    thisoff = [thisoff,switch_off];
                    if switch_on<runtime
                        thison = [thison,switch_on];
                    else
                        thison = [thison,runtime];
                    end
                else
                    thisoff = [thisoff,runtime];
                end
                currenttime = switch_on;
            end
        end
        thison(thison<0) = [];
        thisoff(thisoff<0) = [];
        thison(thison==runtime) = [];
        thisoff(thisoff==runtime) = [];
        if isempty(thison)
            %Do nothing
        elseif isempty(thisoff)
            AdjMat(i,j) = 1;
            AdjMat(j,i) = 1;
        elseif thison(end)>thisoff(end)
            AdjMat(i,j) = 1;
            AdjMat(j,i) = 1;
        else
            %Do nothing
        end
    end
end
end