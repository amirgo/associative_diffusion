rounds = 1000000; % simulation length
N = 500; % population size
K = 6; % number of practices
outsamples = 100; % number of output iterations to report
referenceB = 100; % references for the gap statistic
gsize = 20; % number of gap statistic calculations per simulation (computationally costly)
s = 6; % number of simulations

% initialize output data structures
gaps = zeros(s,gsize);

results1 = zeros(s,outsamples+1);
results2 = zeros(s,outsamples+1);
results3 = zeros(s,outsamples+1);

decays = zeros(s,1);

% run a prallel pool
pobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(pobj)
    pobj = parpool;
end


parfor x=1:s
    fprintf('simulation %d\n',x);
    
    if mod(x,2)==0
        selection = struct('fixed',1,'dist',{'unid',1}); % fully connected
    else
        selection = struct('fixed',0,'network','caveman'); % small world
    end
    decays(x) = 0.2 + 0.6*round(rand); % randomly assign decay to be 0.2 or 0.8

    params = struct('relaxed',0,'groups',5,'decay',decays(x));
    
    [ Rho, Vs ] = GroupSimulateLean(rounds,N,K,params,selection,outsamples);
    results1(x,:) = Rho(:,1)'; % mean abs correlation between Vs
    results2(x,:) = Rho(:,2)'; % mutual information
    results3(x,:) = Rho(:,3)'; % mean interpretative distance
    
    % compute number of culsters using the gap statistic:
    tmp1 = zeros(1,gsize);
    for i=1:gsize
        iix = round(1+i*outsamples/gsize);
       
        v = Vs(:,:,iix);
        v = v./repmat(reshape(max(abs(v')),N,1),1,K);
        try
            eva = evalclusters(v,'kmeans','gap','KList',[1:K*2], 'Distance','Correlation','SearchMethod','firstMaxSE','B',referenceB);
            t = tabulate(eva.OptimalY);
            tmp1(i) = sum(t(:,2)>1); % consider clusters with more than 1 agent
        catch ME            
        end
    end
    
    gaps(x,:) = tmp1;  
end

delete(pobj);