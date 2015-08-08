function g = bn_g(a, b, s)
g = exp(1i*s*b-a*b)./((a-1i*s).*(1+a-1i*s));