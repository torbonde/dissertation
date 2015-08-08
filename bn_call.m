function P = bn_call(model, K)
phi = @(u) terminal_mgf(model, u);
P = exp(-model.r*model.T)*model.S0*bn_p1(K/model.S0, phi);