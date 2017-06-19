function [] = centrality_over_time(input_folder,input_filename)

timestep_for_capture = 900;

timestamp = datestr(now,'yyyymmddTHHMMSS');
structure = '%f %f %f %*s %*s';

iF = ['input/',input_folder];
oF = ['output_',timestamp];

clean_input = strrep(input_filename, '.', '');
dir_ref = [oF,'\',clean_input];
mkdir(dir_ref);

input = [iF,'/',input_filename];

fid = fopen(input);
rawdata = textscan(fid,structure,'Delimiter',',');
fclose(fid);

%==Extract and Clean Data==%
data = cell2mat(rawdata);
data(:,1) = data(:,1)-data(1,1);
lowestID = min(min(data(:,2)),min(data(:,3)));
data(:,2) = data(:,2)-lowestID+1;
data(:,3) = data(:,3)-lowestID+1;
number_rows = size(data,1);
parfor i=1:number_rows
    thisrow = data(i,:);
    col2 = thisrow(1,2);
    col3 = thisrow(1,3);
    if col2 > col3
        thisrow(1,2) = col3;
        thisrow(1,3) = col2;
        data(i,:) = thisrow;
    end
end
all_IDs = [data(:,2); data(:,3)];
all_active = unique(all_IDs);
num_people = size(all_active,1);
data2 = data(:,2);
data3 = data(:,3);
for i=1:num_people
    oldID = all_active(i);
    data2(data2==oldID) = -i;
    data3(data3==oldID) = -i;
end
data(:,2) = -data2;
data(:,3) = -data3;

maxtime = max(data(:,1));
capturepoints = 0:timestep_for_capture:maxtime;

for i=1:length(capturepoints)
    currentcapturetime = capturepoints(i);
    currentconnections = data(data(:,1)==currentcapturetime,2:3);
    
    AdjMat = zeros(num_people);
    
    if size(currentconnections,1)>0
        for j=1:size(currentconnections,1)
            x = currentconnections(j,1);
            y = currentconnections(j,2);
            AdjMat(x,y) = 1;
            AdjMat(y,x) = 1;
        end
    end
    
    connections = sum(AdjMat);
    histimg = figure();
    histogram(connections);
    xlabel('Node Degree');
    ylabel('Frequency');
    titletext = ['Node Degree Frequency - t=',num2str(currentcapturetime)];
    title(titletext);
    imagefilename = [dir_ref,'/','nodedegreehist-t',num2str(currentcapturetime),'.png'];
    print(imagefilename,'-dpng')
    close(histimg);
    
    G = graph(AdjMat);
    graphimg = figure();
    plot(G);
    axis off
    titletext = ['Network Map - t=',num2str(currentcapturetime)];
    title(titletext);
    imagefilename = [dir_ref,'/','networkmap-t',num2str(currentcapturetime),'.png'];
    print(imagefilename,'-dpng')
    close(graphimg);
    
    CCS = figure();
    p = plot(G);
    ucc = centrality(G,'closeness');
    p.NodeCData = ucc;
    colormap jet
    colorbar
    axis off
    titletext = ['Closeness Centrality Scores - t=',num2str(currentcapturetime)];
    title(titletext);
    imagefilename = [dir_ref,'/','closenesscentralityscores-t',num2str(currentcapturetime),'.png'];
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
    titletext = ['Betweenness Centrality Scores - t=',num2str(currentcapturetime)];
    title(titletext);
    imagefilename = [dir_ref,'/','betweenesscentralityscores-t',num2str(currentcapturetime),'.png'];
    print(imagefilename,'-dpng')
    close(BCS);
    
    ECS = figure();
    p = plot(G);
    ucc = centrality(G,'eigenvector');
    p.NodeCData = ucc;
    colormap jet
    colorbar
    axis off
    titletext = ['Eigenvector Centrality Scores - t=',num2str(currentcapturetime)];
    title(titletext);
    imagefilename = [dir_ref,'/','eigenvectorcentralityscores-t',num2str(currentcapturetime),'.png'];
    print(imagefilename,'-dpng')
    close(ECS);
    
    DCS = figure();
    p = plot(G);
    ucc = centrality(G,'degree');
    p.NodeCData = ucc;
    colormap jet
    colorbar
    axis off
    titletext = ['Degree Centrality Scores - t=',num2str(currentcapturetime)];
    title(titletext);
    imagefilename = [dir_ref,'/','degreecentralityscores-t',num2str(currentcapturetime),'.png'];
    print(imagefilename,'-dpng')
    close(DCS);
    
    PCS = figure();
    p = plot(G);
    ucc = centrality(G,'pagerank');
    p.NodeCData = ucc;
    colormap jet
    colorbar
    axis off
    titletext = ['Pagerank Centrality Scores - t=',num2str(currentcapturetime)];
    title(titletext);
    imagefilename = [dir_ref,'/','pagerankcentralityscores-t',num2str(currentcapturetime),'.png'];
    print(imagefilename,'-dpng')
    close(PCS);
end