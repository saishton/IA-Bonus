function [datamatrix] = model4fullopt(nodes,runtime,change)

cut = 20;

preruntime = 20*ones(nodes);
switchon = exprnd(7384.5,nodes);
startthings = switchon-preruntime;

initial = zeros(nodes);

ex1_mu = change;
ln1_mu = 1.8887;
ln2_mu = 0.2536;

%LNpara1 = lognrnd(ln1_mu,sigma_for_mu_and_mean(6.6106,ln1_mu),nodes);
%LNpara2 = lognrnd(ln2_mu,sigma_for_mu_and_mean(1.2887,ln2_mu),nodes);
EXpara1 = lognrnd(ex1_mu,sigma_for_mu_and_mean(30.552,ex1_mu),nodes);
LNpara1 = 6.6106*ones(nodes);
LNpara2 = 1.2887*ones(nodes);
%EXpara1 = 30.552*ones(nodes);

ontimes = struct();
offtimes = struct();

for i=1:nodes-1
    for j=i+1:nodes
        init = initial(i,j);
        currenttime = startthings(i,j);
        if init == 0
            thisoff = [startthings(i,j)];
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
            thison = [startthings(i,j)];
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
        firstonIDX = find(thison>0,1);
        firstoffIDX = find(thisoff>0,1);
        firston = thison(firstonIDX);
        firstoff = thisoff(firstoffIDX);
        thison(thison<0) = [];
        thisoff(thisoff<0) = [];
        thison(thison==runtime) = [];
        if firston > firstoff
            thison = [switchon(i,j),thison];
        end
        ID_ref = sprintf('n%d_n%d', i,j);
        ontimes.(ID_ref) = thison;
        offtimes.(ID_ref) = thisoff;
    end
end
datamatrix = sampleCSV4opt(ontimes,offtimes,nodes,runtime,cut);
end