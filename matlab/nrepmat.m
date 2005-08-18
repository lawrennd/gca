function B = nrepmat(A,dim,number)

% NREPMAT Replicate and tile an array.

% GCA

evaltimes=abs(',A')';
evaltimes=evaltimes*ones(1, number);
evaltimes=char(evaltimes(:)');

eval(['B=cat(dim ' evaltimes ');'])
