function c = log_pde_ko_call(model, option, M, N)
% LOG_PDE_KO_CALL prices a knock-out call option using finite differences on
% the log transformed Black-Scholes PDE.
% 
%   Required are:
%       model : struct with fields
%           sigma : Volatility
%           r     : Risk-free interest
%           S0    : Initial value
%       option : struct with fields
%           T : Time to maturity
%           K : Strike price
%           L : Lower barrier
%           U : Upper barrier
%       M : Number of pieces the domain will be divided into.
%       N : Number of pieces [0,T] will be divided into.
opt.T = option.T;
opt.h = @(x) max(exp(x)-option.K,0);
opt.L = log(option.L);
opt.U = log(option.U);

model.x0 = log(model.S0);
[c, ~] = log_pde_solve(model, opt, M, N);