function comp_val = compare4opt(real_data,gen_data)

numcomps = 2;

comp_val = Inf(numcomps,1);

[~,~,ks2stat] = kstest2(real_data,gen_data);

Z = sort([real_data,gen_data]);
X = sort(gen_data);

idx = zeros(1,length(Z));
parfor i=1:length(Z)
    [~,raw] = find(X<=Z(i),1,'last');
    if isempty(raw)
        idx(i)=0;
    else
        idx(i)=raw;
    end
end

N = length(real_data)+length(gen_data);
sumterm = zeros(1,N-1);
parfor i=1:N-1
    sumterm(i) = ((idx(i)*N-length(real_data)*i)^2)/(i*(N-i));
end
adsum = sum(sumterm);
ad2stat = adsum/(length(real_data)*length(gen_data));

comp_val(1) = ks2stat;
comp_val(2) = ad2stat;

    

