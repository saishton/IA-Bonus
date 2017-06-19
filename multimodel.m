function [] = multimodel(nodes,runtime,number)

timestamp = datestr(now,'yyyymmddTHHMMSS');
dir_ref = ['output_',timestamp];
mkdir(dir_ref);

for runnum=1:number
    
    cut = 20;
    
    lambda_set = 6278.0;
    mu1_set = 3.5348;
    sig1_set = 0.2807;
    mu2_set = 6.3512;
    sig2_set = 1.3688;
    
    preruntime = 20*ones(nodes);
    switchon = exprnd(lambda_set,nodes);
    startthings = switchon-preruntime;
    
    initial = zeros(nodes);
    
   
    EXpara1 = lognrnd(mu1_set,sig1_set,nodes);
    LNpara1 = mu2_set*ones(nodes);
    LNpara2 = sig2_set*ones(nodes);
    
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
    
    filename = sprintf('gen_data_%d.csv', runnum);
    multisampleCSV(ontimes,offtimes,nodes,runtime,cut,dir_ref,filename);
end

end