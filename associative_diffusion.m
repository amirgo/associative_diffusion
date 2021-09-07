function [ Rho, Vs ] = associative_diffusion(rounds,N,K,params,selection,outputSamples)
%   Input: 
%   rounds: number of simulation rounds
%   N: number of agents
%   K: number of behaviors
%   params: a struct with model parameters:
%           relaxed: a binary indicator of relaxed constraint (see online supplement, Part B)
%           groups: numbe of groups, relevnt only if a caveman network
%           decay: the retention parameter (lambda), ranging form 0 to 1
%   selection is a struct that defines how the actor is selected:
%           fixed: binary if based on a fixed dist (set this to 1 if you
%                  don't want social structure to be taken into account)
%           dist: distribution
%           prominence: an N sized vector of weights of prominence
%           network: caveman (specify their groups otherwise no effect) or scalefree
%   outSamples: how many samples to output
% 
%   Output:
%   Rho: statistics by iteration:
%        column 1: preference congruence
%        column 2: mutual information
%        column 3: interpretative distance 
%   Vs: agents' preference vectors by iteration


% initialize parameters to defaults
subgroups = 1;
decay = 0.95;
network_beta = -1;

if nargin>=5
    if isfield(params,'relaxed')
        relaxed = params.relaxed;
    end
    if isfield(params,'groups')
        subgroups = params.groups;
    end
    if isfield(params,'decay')
        decay = params.decay;
    end
end

% calculate interval of output 
if outputSamples
    increment = floor(rounds/outputSamples);
    if (outputSamples*increment)==rounds
        outputSamples = outputSamples + 1;
    else
        outputSamples = outputSamples + 2;
    end
else
    increment = 1;
    outputSamples = rounds;
end

% divide into groups
groups = [];
for i=1:subgroups
    groups = [groups,i*ones(1,round(N/subgroups))];
end
if length(groups)<N
    groups = [groups,i*ones(1,N-length(groups))];
else
    groups = groups(1:N);
end

% G = unique(groups);

S = zeros(N,N); % The social network
R = ones(K,K,N); % Agents' relational matrices


prominence = ones(1,N); % The prominence of each agent, initialized to 1
fixedSelection = 1;
% conformity = -1*ones(1,N);

% determine selection mode:
if nargin>=5
    if isfield(selection,'fixed')
        fixedSelection = selection.fixed;
    end
    if isfield(selection,'prominence')
        prominence = selection.prominence;
    elseif isfield(selection,'dist')
        prominence = generate_dist({selection.dist},1,N);
    end
    if isfield(selection,'network')
        S = generate_network({selection.network},groups,network_beta);
    end
end


% set self-associations to 0
for i=1:N
    R(:,:,i) = ones(K)-eye(K);
end

Rho = zeros(outputSamples,3);
Vs = zeros(N,K,outputSamples);

V = generate_dist({'unif',-1,1},N,K); % initialize preferences to random with uniform probability

indx = triu(ones(N),1)==1;

six = 0;
for x=1:rounds
    if ~mod(x,1000)
        fprintf('.'); % show progress
    end
    
    % choose two random agents
    observer = ceil(N*rand(1,1)); % choose observer at random
    if ~fixedSelection % if selection isn't fixed, set selection probability to the agent's network ties
        prominence = S(observer,:);
    end
    
    % get the actor as a function of the observer, with probability porportional to prominence
    actor = get_actor(observer,prominence);

    [i,j] = create_behavior_lean(V(actor,:)); % actor enacts 2 practices

    % update observer's associative matrix and preference, conditional on constraint satisfaction:
    [V(observer,:),R(:,:,observer)] = observe_behavior(i,j,R(:,:,observer),V(observer,:),relaxed);
    
    R(:,:,observer) = decay*R(:,:,observer); % apply decay

    % compute statistics per output interval:
    if x==1 || mod(x,increment)==0 || x==rounds
        six = six+1;
        

        Ccur = corr(V');
        Ccur = Ccur.*(1-eye(N));
        if all(all(isnan(Ccur)))
            Ccur(isnan(Ccur)) = 1;
        end
        
        Rho(six,1) = mean(abs(Ccur(indx)));
        [mi_] = mutual_information_lean(V);
        Rho(six,2) = mi_;
        [F_, rho_] = forbenius(R);
        Rho(six,3) = rho_;
        

        
        Vs(:,:,six) = V;
    
    end
    

end
fprintf('\n');
end


function [ i,j ] = create_behavior_lean(Atts)
%   enact to practices on the basis of preferences


v = Atts;
p = get_behavior_probabilities_lean(v);


if sum(isnan(p)==1)>0
    p = 0;
end
i = discreteinvrnd(p,1,1);

% now choose j that is not i
p(i) = 0;
if all(p==0)
    p = p+1;
    p(i) = 0;
end
p = p./sum(p);


j = discreteinvrnd(p,1,1);

while j==i
    j = discreteinvrnd(p,1,1);
end

end


function [ p ] = get_behavior_probabilities_lean(V)
%   exponentiate and normalize probabilities to 1

K = size(V,2);
p = exp(V);

p(p==Inf) = realmax;
if all(p==0)
    p(p==0) = 1;
end
p = p./repmat(sum(p,2),1,K);


end




function [ newAtts,newR ] = observe_behavior(i,j,R,Atts,relaxed)
%   update R and V once i and j are observed

newR = R;
newR(i,j) = R(i,j)+1;
newR(j,i) = R(j,i)+1;

newAtts = update_attitudes(Atts,R,newR);

% if constraint satisfaction isn't changed, keep old preferences
if relational_compare(newR,newAtts,relaxed) <= relational_compare(newR,Atts,relaxed)
    newAtts = Atts;
end

newAtts(newAtts==Inf) = realmax;
newAtts(newAtts==-Inf) = -realmax;

end



function [ rho, pearson, frobenius ] = relational_compare(R,V,relaxed)
%   compute constraint satisfaction

n = size(R,1);

dV = abs(repmat(V',1,n)-repmat(V,n,1));

dV = maxminnorm(dV);

r = maxminnorm(R);


mP = ones(n,n);
if relaxed
    p = reshape(V>0,1,n);
    mP = max(repmat(p,n,1),repmat(p',1,n));
end
D = mP.*abs(r-dV);

indx = triu(ones(n),1)==1;
rho = sum(D(indx));
rho = rho*2/((n-1)*n);

if nargout>1
    pearson = corr(r(indx),dV(indx));
end

if nargout>2
    frobenius = sqrt(sum(diag((r-dV)*((r-dV)'))));
end

end


function [R] = maxminnorm(X)
%   standardize R


mn = min(min(X));
mx = max(max(X));

if mn~=mx
    R = X./mx;
else
    R = X-X;
end
end



function [ newAtts ] = update_attitudes(Atts,R0,R1)
%   update the weaker preference randomly

newAtts = Atts;


v = Atts;

n = size(R0,1);

if size(R0,1)~=size(R1,1)
    fprintf('R0 and R1 should have the same size\n');
    return;
end

indx = triu(ones(n),1)==1;

r0 = R0(indx);
r1 = R1(indx);

v = (v-min(v));
v = v./max(v);

v_hat = mean(v);


diff = find(r0 ~= r1);

v_new = Atts;

for ix=1:length(diff)
    [i,j] = getSubscript(diff(ix));
    if abs(v(i)-v_hat) >= abs(v(j)-v_hat)
        x = j;
    else
        x = i;
    end
    
    v_new(x) = v_new(x)+randn;
    
end

newAtts = v_new; 

end

function [i,j] = getSubscript(indx)

y = floor(0.5*(sqrt(1+8*indx)-1));
i = 0.5*y*(y+1);

if mod(indx,i) || (i==1 && indx==2)
    j = y+2;
    if i>1
        i = mod(indx,i);
    end
else
    j = y+1;
    i = y;   
end

end


function [MI] = mutual_information_lean(V)
% mutual information

N = size(V,1);
K = size(V,2);

p = get_behavior_probabilities_lean(V);

J= zeros(K);
p_1 = zeros(K,1);
p_2 = zeros(1,K);
for i=1:N
    j = condprob(p(i,:)');
    J = J + j./N;
    p_1 =  p_1 + sum(j,2)./N;
    p_2 =  p_2 + sum(j,1)./N;
end

P = p_1*p_2;


indx = eye(K)==0;


MI = J.*log2(J./P);

MI = sum(MI(indx));

end


function [ actor ] = get_actor(observer,probabilities)
%   get an actor that isn't observer


probabilities(probabilities<0) = 0;
p = probabilities;
p(observer) = 0;
if(all(p==0))
    p = p+1;
end
actor = discreteinvrnd(p,1,1);

end


function condP = condprob(p)

K = length(p);
condP = zeros(K);

for i=1:K
    ix = [1:K]~=i;
    px = p(i);
    py = p(ix);
    
    py = py./repmat(sum(py),K-1,1);
    condP(i,ix) = (repmat(px,K-1,1).*py)';
end

end











