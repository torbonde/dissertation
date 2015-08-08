function phi = terminal_mgf(model, u)
phi = exp((model.r - 1/2*model.sigma^2)*u + 1/2*model.T*model.sigma^2*u.^2);