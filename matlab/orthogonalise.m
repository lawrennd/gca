function W = orthogonalise(W);

% ORTHOGONALISE Orthogonalise a principal subspace.

% GCA

[U, V] = eig(W*W');
W = U'*W;
