function [ POST_PROB, PEPTIDE_IDX, NUM_PROTEOFORMS ] = Generate_Posterior( COUNTS, PI_PROBS, P_CONFIGS, COUNTS_IDX )
%This function is used to generate posterior probability for BP-Quant

    % INPUTS:
        % COUNTS = a vector of size (Nu x 1) that represents the count of 
            % each unique signature (Nu total) in PROTEIN_SIG
        % PI_PROBS = a vector of size (Nu x 1) that represents the
            % background frequency expected for each signature (Nu total)
        % P_CONFIGS = a vector of size (Nk x Nu) that represents the
            % configuration across the Nu signatures where a 0 indicates
            % there is no proteoform and 1 indicates that there is one.
            
    % REQUIRED OUTPUTS:
        % POST_PROB = a vector of size (Nk x Nu) that represents the
            % probability for each of the Nk configurations
            
    % Bobbie-Jo Webb-Robertson (10/22/13)


%% DEFINE SIZES
Nu = length(COUNTS);
N_PEPS = sum(COUNTS);
N_CONFIGS = size(P_CONFIGS,1);

%% DEFINE CPM
X = zeros(Nu,2);
for i = 1:Nu
    X(i,2) = binocdf(COUNTS(i)-1,N_PEPS,PI_PROBS(i));
end
X(:,1) = 1-X(:,2);
    
%% GENERATE POSTERIOR
POST_PROB = ones(N_CONFIGS,1);
for i = 1:N_CONFIGS
    x = P_CONFIGS(i,:);
    for j = 1:Nu
        prior = PI_PROBS(j)^x(j)*(1-PI_PROBS(j))^(1-x(j));
        POST_PROB(i) = POST_PROB(i)*prior*X(j,x(j)+1);
    end
end
POST_PROB = POST_PROB ./ sum(POST_PROB);

%% DETERMINE THE NUMBER OF PROTEOFORMS AND MAP THEM TO THE APPOPRIATE PEPTIDE ORDER IN PROTEIN_SIG
[~,b] = max(POST_PROB);
P_tmp = P_CONFIGS(b,:);
if sum(P_tmp) == 0
    PEPTIDE_IDX = ones(N_PEPS,1);
    NUM_PROTEOFORMS = 1;
else
    PEPTIDE_IDX = zeros(N_PEPS,1);
    NUM_PROTEOFORMS = sum(P_tmp);
    t = find(P_tmp > 0);
    k = 0;
    for i = 1:NUM_PROTEOFORMS
        k = k+1;
        PEPTIDE_IDX(COUNTS_IDX==t(i)) = i;
    end
end

end

