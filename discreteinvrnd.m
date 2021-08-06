function X = discreteinvrnd(p,m,n)

% MATLAB's discreteinvrnd function

    X = zeros(m,n);
    for i = 1:m*n
        u = rand;
        I = find(u < cumsum(p));
        X(i) = min(I);
    end
end
