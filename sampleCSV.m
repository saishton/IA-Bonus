function [] = sampleCSV(ontimes,offtimes,nodes,runtime,sampletime)

timestamp = datestr(now,'yyyymmddTHHMMSS');
dir_ref = ['output_',timestamp];
mkdir(dir_ref);

numint = floor(runtime/sampletime)+1;
timesteps = 0:sampletime:runtime;
massivematrix = [];

for i=1:nodes-1
    for j=i+1:nodes
        ID_ref = sprintf('n%d_n%d', i,j);
        thison = ontimes.(ID_ref);
        thisoff = offtimes.(ID_ref);
        
        thisindicator = zeros(1,numint);
        parfor i=1:numint
            currenttime = (i-1)*sampletime;
            if sum(thison<=currenttime & thisoff>=currenttime)
                thisindicator(i) = 1;
            end
        end
        thistimes = thisindicator.*timesteps;
        thistimes(thistimes==0) = [];
        thistimes = thistimes';
        n = length(thistimes);
        thisc2 = i*ones(n,1);
        thisc3 = j*ones(n,1);
        thisc4 = ones(n,1);
        thisc5 = ones(n,1);
        thisblock = [thistimes,thisc2,thisc3,thisc4,thisc5];
        massivematrix = [massivematrix;thisblock];
    end
end
[~,idx] = sort(massivematrix(:,1));
sortedbytime = massivematrix(idx,:);

filename = 'generated_data.csv';
filepath = [dir_ref,'/',filename];
csvwrite(filepath,sortedbytime)