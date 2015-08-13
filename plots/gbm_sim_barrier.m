rng('default');
close all; clear; clc;
M = 200;
N = 100;

L = [45, 40, 35, 30];
U = [55, 60, 65, 70];

Bmin = min(L)-5;
Bmax = max(U)+5;

filename = @(i) ['../../figures/data/gbm_sim_barrier_' num2str(i) '.csv'];

for i = 1:4
   subplot(2,2,i);
   [ts, X] = gbm_sim(M, N);
   Y = nan_hits(X, L(i), U(i));
   hold on
   plot(ts,Y);
   plot([0 1], [L(i) L(i)], ':k')
   plot([0 1], [U(i) U(i)], ':k')
   hold off
   ylim([Bmin, Bmax]);
   xlim([0 1]);
   legend(['Barrier hits: ' num2str(sum(isnan(Y(:,end)))) ' / ' num2str(M)]);
   
   %csvwrite(filename(i), [ts', Y']);
end