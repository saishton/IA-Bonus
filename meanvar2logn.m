function [mu,sig] = meanvar2logn(m,v)

mu = log((m^2)/sqrt(v+(m^2)));
sig = sqrt(log((v/(m^2))+1));

end