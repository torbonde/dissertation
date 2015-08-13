function [ts, Y] = gbm_sim(M,N)
S0 = 50;
r = 0.05;
sigma = 0.2;
T = 1;

dt = T/N;

Y = zeros(M,N+1);
Y(:,1) = S0;
for i = 2:N+1
    Y(:,i) = Y(:,i-1).*exp((r-sigma^2/2)*dt + sigma*sqrt(dt)*randn(M,1));
end
ts = 0:dt:T;