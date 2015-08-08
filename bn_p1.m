function P1 = bn_p1(K, phi)
a=9; % optional?
g = @(s) bn_g(a, log(K), s);
P1 = real(1/(2*pi)*integral(@(s) phi(1+a-1i*s).*g(s), -Inf, Inf)); % `real` just to satisfy Matlab