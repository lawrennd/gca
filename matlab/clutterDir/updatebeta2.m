function beta = updatebeta2(sbar, Sigma_s, XU, expD, expD2, V);
A = expD*V';
a = 0.5;
b = 0.5;
ndata = size(XU, 1);
dataDim = size(XU, 2);
beta = 0;
% Need to use expected values for A!!!!

for i = 1:dataDim
  for n = 1:ndata
    %     beta = beta + XU(n, i)*XU(n, i) ...
    % 	   - 2*XU(n, i)*A(i, :)*sbar(n, :)' ...
    % 	   + (A(i, :)*sbar(n, :)')*(sbar(n, :)*A(i, :)') ...
    % 	   + A(i, :)*Sigma_s(:, :, n)*A(i, :)';
    beta = beta + (XU(n, i) - A(i, :)*sbar(n, :)')^2 + A(i, :)*Sigma_s(:, :, n)*A(i, :)';
    
  end
end
beta = (a + 0.5*dataDim*ndata)/(b+0.5*beta);

