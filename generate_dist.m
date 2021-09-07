function D = generate_dist(dist,N,K,prefix)

if nargin<4
    prefix = '';
end

if strcmp(dist{1},'ratio')
    ratio = dist{2};
    strength = dist{3};
    ix = randperm(N*K);
    ix = ix(1:round(N*K*ratio));
    D = strength*ones(N,K);
    D(ix) = 1-strength;
else
    disttext = sprintf('''%s''',dist{1});
    for i=2:length(dist)
        disttext = sprintf('%s,%.10f',disttext,dist{i});
    end

    evaltext = sprintf('D = %srandom(%s,N,K);',prefix,disttext);
    eval(evaltext);

    if strcmp(dist{1},'wbl')
        D = D-min(D);
        D = D./max(D);
    end
end

end