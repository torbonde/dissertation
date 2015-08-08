function P = bn_uo_call(model, K)
phi = @(u) terminal_mgf(model, u);
psi = @(u, v) terminal_max_mgf(model, u, v);
B = log(model.U/model.S0);
P = exp(-model.r*model.T)*model.S0*bn_p2(K, B, phi, psi);