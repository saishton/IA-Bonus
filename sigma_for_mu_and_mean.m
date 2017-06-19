function [sigma] = sigma_for_mu_and_mean(m,mu)

if log(m)>mu
    sigma = sqrt(2*(log(m)-mu));
else
    sigma = 0;
end
end