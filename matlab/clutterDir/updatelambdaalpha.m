function alpha = updatelambdaalpha(lambda, sigma2_lambda);

a = 1e-3;
b = 1e-3;

dataDim = length(lambda);
alpha = zeros(1, dataDim);
for i = 1:dataDim
  alpha(i) = (0.5+a)/(b+0.5*(lambda(i)*lambda(i) + sigma2_lambda(i)));
end


  




