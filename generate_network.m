function S = generate_network(type,groups,beta)

S = zeros(length(groups));
if strcmp(type,'caveman')
    g = unique(groups);
    for i=1:length(g)
        ix = groups==g(i);
        S(ix,ix) = 1;
    end
    S = S-eye(size(S));
    prop = 0.1*length(find(S==1))/2;
    S=sym_generate_srand(S,round(prop));
elseif strcmp(type,'full')
    S = ones(length(groups));
    S = S-eye(size(S));
end
