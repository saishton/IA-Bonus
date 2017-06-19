function [transprobs] = transitionProb(init,runtime)

sims = 1E3;

simssol = zeros(1,sims);

parfor i=1:sims
    am = modelAM(init,runtime);
    AL = sum(sum(am))/2;
    simssol (i) = AL;
end

uni = unique(simssol);
count = zeros(1,length(uni));

parfor i=1:length(uni)
    count(i) = sum(simssol==uni(i));
end
probs = count/sum(count);
transprobs = [uni;probs];
