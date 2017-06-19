function [AdjMat] = modelAMfigs(initMat,runtime)

nodes = size(initMat,1);

timestamp = datestr(now,'yyyymmddTHHMMSS');
dir_ref = ['output_',timestamp];
mkdir(dir_ref);

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
triangles = trace(AdjMat^3)/6;
triples = (sum(sum(AdjMat^2))-trace(AdjMat^2))/2;
GCC = 3*triangles/triples;

connections = sum(AdjMat);
histimg = figure();
histogram(connections);
xlabel('Node Degree');
ylabel('Frequency');
imagefilename = [dir_ref,'/','nodedegreehist.png'];
print(imagefilename,'-dpng')
close(histimg);

G = graph(AdjMat);
graphimg = figure();
plot(G);
axis off
title('Network Map');
imagefilename = [dir_ref,'/','networkmap.png'];
print(imagefilename,'-dpng')
close(graphimg);

CCS = figure();
p = plot(G);
ucc = centrality(G,'closeness');
p.NodeCData = ucc;
colormap jet
colorbar
axis off
title('Closeness Centrality Scores');
imagefilename = [dir_ref,'/','closenesscentralityscores.png'];
print(imagefilename,'-dpng')
close(CCS);

BCS = figure();
p = plot(G);
wbc = centrality(G,'betweenness');
n = numnodes(G);
p.NodeCData = 2*wbc./((n-2)*(n-1));
colormap(flip(autumn,1));
colorbar
axis off
title('Betweenness Centrality Scores');
imagefilename = [dir_ref,'/','betweenesscentralityscores.png'];
print(imagefilename,'-dpng')
close(BCS);

ECS = figure();
p = plot(G);
ucc = centrality(G,'eigenvector');
p.NodeCData = ucc;
colormap jet
colorbar
axis off
title('Eigenvector Centrality Scores');
imagefilename = [dir_ref,'/','eigenvectorcentralityscores.png'];
print(imagefilename,'-dpng')
close(ECS);

DCS = figure();
p = plot(G);
ucc = centrality(G,'degree');
p.NodeCData = ucc;
colormap jet
colorbar
axis off
title('Degree Centrality Scores');
imagefilename = [dir_ref,'/','degreecentralityscores.png'];
print(imagefilename,'-dpng')
close(DCS);

PCS = figure();
p = plot(G);
ucc = centrality(G,'pagerank');
p.NodeCData = ucc;
colormap jet
colorbar
axis off
title('Pagerank Centrality Scores');
imagefilename = [dir_ref,'/','pagerankcentralityscores.png'];
print(imagefilename,'-dpng')
close(PCS);
end