function P = smc(model, N, M)
% SMC Runs the sequential Monte Carlo algorithm.
%
%   P = SMC(model, N, M) estimates the price P of a barrier option at time
%   0 using the sequential Monte Carlo estimator under the Black-Scholes
%   model with parameters specified  by the struct model. M trajectories
%   are simulated at N intermediate time-steps to estimate the time 0-price P.
%   The struct model must contain the fields
%       T      : Time to maturity
%       S0     : Spot price
%       r      : Risk-free interest
%       sigma  : Volatility
%       payoff : Function of terminal value, describing payoff at maturity.
%   and at least one of
%       L      : Lower boundary
%       U      : Upper boundary

P = 0;
dt = model.T/N;
dW = randn(M,N);
G = zeros(M,N);
E = zeros(N,1);
Snm1 = repmat(model.S0, M, 1);

% Check boundaries
if isfield(model, 'L') && isfield(model, 'U')
    ind = @(s) model.L < s & s < model.U;
    
    cutoff = 3;
    f = @(x,y) 1 - double_barrier_exit_prob(x, y, model.L, model.U, model.sigma, dt, cutoff);
elseif isfield(model, 'L')
    ind = @(s) model.L < s;
    f = @(x,y) 1 - single_barrier_exit_prob(x, y, model.L, model.sigma, dt);
elseif isfield(model, 'U')
    ind = @(s) s < model.U;
    f = @(x,y) 1 - single_barrier_exit_prob(x, y, model.U, model.sigma, dt);
else
    error('No boundary specified in model.');
end

for n = 1:N
    Sn = Snm1.*exp((model.r - model.sigma^2/2)*dt + model.sigma*sqrt(dt)*dW(:,n));
    G(:,n) = ind(Sn).*f(Snm1, Sn);
    
    % If all trajectories crossed a barrier, the estimator is null
    if ~any(G(:,n))
        return
    end
    
    % Acceptance/rejection
    U = rand(M,1);
    reject = G(:,n) < U;
    
    % Recycle
    w = G(:,n)/sum(G(:,n));
    num_rejects = sum(reject);
    Sn(reject, :) = sample_weighted(num_rejects, Sn, w);
    
    % Estimate expectation sequentially
    E(n) = sum(G(:,n))/M;
    
    % Store resampled Sn for next step.
    Snm1 = Sn;
end

payoff = sum(model.payoff(Sn))/M;
P = exp(-model.r*model.T)*payoff*prod(E);