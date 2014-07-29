function [ COUNTS, COUNTS_IDX, PI_PROBS, SIGS ] = Parameter_Trends( PROTEIN_SIG, PI_NOT )
%This function is used to generate signature parameters for BP-Quant
%   The function inputs the signature matrix for a protein and generates
%   the unique signatures and their counts, as well as the expected
%   background frequencies.  The data is sorted with the [zero] signature
%   first (if present) and the remaining signatures in order of their
%   counts.  The ordering is simply to make outputs easier to read.

    % INPUTS:
        % PROTEIN_SIG = (Ns x Nt) matrix that contains values of only -1, 0
            % or 1.  Ns represents the number of peptides and Nt is the
            % length of the signature (number of statistical comparisons)
        % PI_NOT = a single value between 0 and 1 that represents the
            % likelihood of observing an all zero trend under the
            % assumption that there is no statistical difference between
            % groups.
            
    % REQUIRED OUTPUTS:
        % COUNTS = a vector of size (Nu x 1) that represents the count of 
            % each unique signature (Nu total) in PROTEIN_SIG
        % COUNTS_IDX = a vector of length sisze (Ns x 1) that identifies
            % which peptides are associated with each signature
        % PI_PROBS = a vector of size (Nu x 1) that represents the
            % background frequency expected for each signature (Nu total)
    % OPTIONAL OUTPUTS:
        % SIGS = a matrix of size (Nu x Nt) that gives the signature in
            % order of the counts.  This is primarily just needed for
            % output that is optional
            
    % Bobbie-Jo Webb-Robertson (10/22/13)
    
%% ERROR CHECKING AND SIZE DEFINITIONS
a = unique(PROTEIN_SIG);
if length(a)>3
    error('The descriptors of the protein signature can only be a -1, 0 or 1')
end
if max(size(PI_NOT))>1 || PI_NOT<0 || PI_NOT>1
    error('The background frequency of the zero signature must be a single number between 0 and 1')
end
Ns = size(PROTEIN_SIG,1);

%% IDENTIFY AND COUNT UNIQUE SIGNATURES
[SIGS,~,t] = unique(PROTEIN_SIG,'rows');
Nu = max(t);
COUNTS = zeros(Nu,1);
for i = 1:Nu
    COUNTS(i) = sum(t==i);
end

%% SORT TRENDS IN DESCENDING ORDER OF COUNTS
[COUNTS,v] = sort(COUNTS,'descend');
SIGS = SIGS(v,:);

%% SET ZERO SIGNATURE FIRST IF PRESENT AND ACCORDINGLY LET THE PROBABILITY VALUES
t = find(sum(abs(SIGS),2)==0);
if isempty(t)
    % no zero signature
    PI_PROBS = repmat((1-PI_NOT)./Nu,Nu,1);
else
    SIGS = [SIGS(t,:);SIGS(1:t-1,:);SIGS(t+1:Nu,:)];
    COUNTS = [COUNTS(t);COUNTS(1:t-1);COUNTS(t+1:Nu)];
    PI_PROBS = repmat((1-PI_NOT)./(Nu-1),Nu,1); 
    PI_PROBS(1) = PI_NOT;
end
COUNTS_IDX = zeros(Ns,1);
for i = 1:Nu    
    COUNTS_IDX((sum(abs(PROTEIN_SIG - repmat(SIGS(i,:),Ns,1)),2)==0)==1) = i;
end

    
