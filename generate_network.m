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
elseif strcmp(type,'scalefree')
    if beta==-1
        beta = 1;
    end
    p = exp(beta*(1:length(groups)));
    p = p./sum(p);
    alpha = 0;
    while alpha<2 || alpha>3
        S = zeros(length(groups));
        for i=1:length(groups)
            p_ = p;
            p_(i) = 0;
            for j=1:6
                q = discreteinvrnd(p_);
                S(i,q) = 1;
                p_(q) = 0;
            end
        end
        [alpha]=powerlaw_fit(sum(S),'finite');
    end
end
