function [F, rho] = forbenius(R)

N = size(R,3);
K = size(R,1);

n = 2/(K*(K-1));

F = 0;
rho = 0;

for i=1:N-1
    for j=i+1:N
        r1 = R(:,:,i);
        r2 = R(:,:,j);
        r1 = r1./max(max(r1));
        r2 = r2./max(max(r2));
        F = F + sqrt(trace((r1-r2)*((r1-r2)')))/K;
        rho  = rho + n*sum(sum(abs(r1-r2)));
    end
end

F = 2*F/((N-1)*N);
rho = 2*rho/((N-1)*N);

end
