function psi = terminal_max_mgf(model, u, v)
normcdf_ = @(x) 0.5*(Faddeeva_erfc(-x/sqrt(2)));
psibar = @(u, v) 2*normcdf_(sqrt(model.T)*(u+v)).*exp(1/2*model.T*(u+v).^2).*(u+v)./(2*u+v) ...
    + normcdf_(1/2*sqrt(model.T)*v).*exp(1/8*model.T*v.^2).*(2*u.*exp(1/4*model.T*u.*(2*u+v))./(2*u+v) - 1);
theta = model.r/model.sigma - model.sigma/2;
psi = exp(-1/2*theta^2*model.T)*psibar(u*model.sigma + theta, v*model.sigma);