function [] = model(nodes,runtime)

preruntime = 1000;

initial = rand(nodes)>0.995;

LNpara1 = lognrnd(1.1664E0,3.3155E-2,nodes);
LNpara2 = lognrnd(-7.2487E-1,2.4700E-1,nodes);
EXpara1 = lognrnd(3.5018E0,1.4879E-1,nodes);

ontimes = struct();
offtimes = struct();

for i=1:nodes-1
    for j=i+1:nodes
        init = initial(i,j);
        currenttime = -preruntime;
        if init == 0
            thisoff = [-preruntime];
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
                    thison = [thisoff,runtime];
                end
                currenttime = switch_off;
            end
        elseif init == 1
            thisoff = [];
            thison = [-preruntime];
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
            thison = [0,thison];
        end
        ID_ref = sprintf('n%d_n%d', i,j);
        ontimes.(ID_ref) = thison;
        offtimes.(ID_ref) = thisoff;
    end
end
sampleCSV(ontimes,offtimes,nodes,runtime,20);
end