rng(60,'twister');
close all; clear; clc;
colors = [1 0 1; ...
          0 1 1; ...
          1 0 0; ...
          0 1 0; ...
          0 0 1; ...
          0 0 0];
[M, ~] = size(colors);
N = 3;
L = 45;
U = 55;

S0 = 50;
r = 0.05;
sigma = 0.2;
T = 1;
dt = T/N;
G = @(s) L < s & s < U;
figure
h = gcf;
hold on
plot([1 N+1], [L L], '--k');
plot([1 N+1], [U U], '--k');

S = zeros(M,N+1);
S(:,1) = S0;
for n = 2:N+1
    S(:,n) = S(:,n-1).*exp((r-sigma^2/2)*dt + sigma*sqrt(dt)*randn(M,1));
    for m = 1:M
        plot(n-1:n, S(m,n-1:n),'-o','Color',colors(m,:));
    end
    wtmp = G(S(:,n));
    w = wtmp/sum(wtmp);
    S(~wtmp,n) = sample_weighted(sum(~wtmp),S(:,n),w);
end

% Get data and save it
% axesObjs = get(h, 'Children');
% dataObjs = get(axesObjs, 'Children');
% xdata = get(dataObjs, 'XData');
% ydata = get(dataObjs, 'YData');
% 
% file = '../../figures/data/resample.csv';
% 
% len = length(xdata)-2;
% Y = zeros(2,2*len);
% Y(:) = nan;
% for l = 1:len
%     y = ydata{l};
%     x = xdata{l};
%     Y(:,l) = x';
%     Y(:,len+l) = y';
% end
% 
% csvwrite(file, Y);