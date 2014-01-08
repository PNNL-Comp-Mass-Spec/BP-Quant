function [ P_CONFIGS ] = Proteoform_Configurations( Nu )
%This function is used to generate signature parameters for BP-Quant
%   The function takes the total count of non-zero signatures and generates
%   the possible proteoform configurations.  The configurations are sorted
%   starting with the [zero] signature.  The ordering is simply to make 
% outputs easier to read.

    % INPUTS:
        % Nu = A single number representing the total number of unique
            % signatures
            
    % REQUIRED OUTPUTS:
        % P_CONFIGS = a vector of size (Nk x Nu) that represents the
            % configuration across the Nu signatures where a 0 indicates
            % there is no proteoform and 1 indicates that there is one.
            
    % Bobbie-Jo Webb-Robertson (10/22/13)
    
%% ERROR CHECKING
if max(size(Nu)) > 1 || Nu-round(Nu) ~= 0
    error('Nu must be a simple integer value representing the number of unique signatures')
end
    
%% BUILD MATRIX OF POSSIBLE PROTEOFORMS
n_combos = 2^Nu;
nr = ones(1,Nu);
for i = 2:Nu
    nr(i) = nr(i-1)*2;
end
k = Nu;
T = [zeros(nr(k),1);ones(nr(k),1)];    
for i = 2:Nu
    k = k-1;
    T(:,i) = repmat([zeros(nr(k),1);ones(nr(k),1)],nr(i),1);
end
    
%% SORT THE CONFIGURATIONS (PRIMARILY JUST FOR EASY TO READ OUTPUT) 
P_CONFIGS = zeros(n_combos,Nu);
k = 1;
for i = 1:Nu
    t = T(sum(T,2)==i,:);
    [t,index] = sortrows(t);
    nt = length(index);
    for j = 0:nt-1
        k = k+1;
        P_CONFIGS(k,:) = t(nt-j,:);
    end
end
    
    
    
end

