function [pdf] = explognconvhist(t)

%Set parameters for exponential and log-normal distributions
mu = 6.3512;
sigma = 1.3688;
lambda = lognrnd(3.5348,0.2807);

Xstep = 1E0;            %PDF slice width
maxn = 1+ceil(t/20);    %Maximum values of n
Xmax = 1E4;             %Last slice for PDFs
N = 1E4;                %Number of MC simulations
M = 1E6;                %Number of each variable in each simulation

X = 0:Xstep:Xmax;
tidx = find(X==t);

pdf = zeros(1,maxn+1);

lnpdf = lognpdf(X,mu,sigma);
lnpdf = lnpdf/sum(lnpdf);       %Normalise PDF
expdf = exppdf(X,lambda);
expdf = expdf/sum(expdf);       %Normalise PDF
copdf = conv(lnpdf,expdf);      %Create PDF for EX+LN

l = length(copdf);
CCDF = zeros(1,l);
parfor i=1:l
    CCDF(i) = 1-sum(copdf(1:i));    %Create CCDF needed
end

thispdf = copdf;
for i=2:maxn+1
    pdf(i) = sum(thispdf(1:tidx).*fliplr(CCDF(1:tidx))); %Convolved*CCDF
    thispdf = conv(thispdf,copdf);                       %Next convolve
end
pdf(1)=1-sum(pdf);      %Calculate probability n=0, by inverse

num=zeros(1,N);
%Begin MC simulation
parfor k=1:N
dt=exprnd(lambda,1,M)+lognrnd(mu,sigma,1,M);
time=cumsum(dt);
index=find(time>t,1);
num(k)=min(index)-1;
end

freq = zeros(1,maxn+1);
parfor j=1:max(num)
    freq(j)=length(find(num==j-1))/N;
end

%Plot on graph
mini = 1E-4;
n=0:maxn;
histimage = figure();
hold on
bar(n,pdf,'hist');
plot(n,freq,'x-');
xlabel('n');
text = ['Prob[N(',num2str(t),')=n]'];
ylabel(text);
histmax = find(pdf>mini,1,'last');
freqmax = find(freq>mini,1,'last');
lastidx = max(histmax,freqmax)+5;
xlim([-0.5 lastidx+0.5])
titleText = ['Parameters: \lambda=',num2str(lambda),', \mu=',num2str(mu),', \sigma=',num2str(sigma)];
title(titleText);
legend('Discrete PDF using Convolve','Monte-Carlo Simulation');
set(gca,'FontSize',18);

end
