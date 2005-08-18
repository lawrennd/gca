function beta = updatebeta(model, X)

% UPDATEBETA Update the noise variance in the variational algorithm.

% GCA

expectedOutput = (X - model.sBar*model.A');

if model.FANoise % Factor analysis type noise model
  beta = zeros(1, model.dataDim);
  for i = 1:model.dataDim
    beta(i) = beta(i) + expectedOutput(:, i)'*expectedOutput(:, i);
    for n = 1:model.numData
      beta(i) = beta(i) + model.A(i, :)*model.Sigma_s(:, :, n)*model.A(i, :)';
    end
    beta(i) = model.numData/beta(i);
  end
else
  beta = 0;
  
  for i = 1:model.dataDim
    beta = beta + expectedOutput(:, i)'*expectedOutput(:, i);
    for n = 1:model.numData
      beta = beta + model.A(i, :)*model.Sigma_s(:, :, n)*model.A(i, :)';
    end
  end
  beta = model.numData*model.dataDim/beta;
end