function sigma2_tau = update_sigma2_tau(tau, nu, a, b);


ndata = size(tau, 1);


sigma2_tau = (a + nu/2*ndata)./(sum(tau, 1).*nu/2 + b);
