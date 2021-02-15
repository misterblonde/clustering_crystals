function [K] = pmkernel(Graphs, L, d)
% Computes the pyramid match graph kernel for a set of unlabeled graphs
%
% Input:
%   Graphs: a 1 x N array of graphs
%   L: number of levels
%   d: dimensionality of node embeddings
%
% Output:
%   K: the N x N kernel matrix
%   runtime: runtime in seconds
%
% Copyright (c) 2016, Giannis Nikolentzos
% nikolentzos@aueb.gr

if nargin < 3
    d = 6;
end

if nargin < 2
    L = 4;
end

N = size(Graphs,1);
%tic;MUTAG

% Compute embeddings
disp('Computing embeddings...');
Us = cell(1,N);
for i=1:N
    n = size(Graphs(i).am, 1);
    [U, ~] = eigs(Graphs(i).am, min(n, d));
    U = abs(U);
    Us{i} = U;
end

% Histogram creation
disp('Creating histograms...');

Hs = cell(N,L);
fprintf(1,'%s','Progress: 0%');
for i=1:N
    last_percent = floor(100*(i-1)/N);
    cur_percent = floor(100*i/N);
    
    if last_percent ~= cur_percent
        last_str = sprintf('%d%s',last_percent,'%');
        fprintf(1,repmat('\b', 1, length(last_str)));
        fprintf(1,'%d%s',cur_percent,'%');
    end
    
    for j=1:L
        l = 2^(j-1);
        D = zeros(d, l);
        T = ceil(Us{i}*l);
        T(T==0) = 1;
        for p=1:size(Us{i},1)
            for q=1:size(Us{i},2)
                D(q,T(p,q)) = D(q,T(p,q)) + 1;
            end
        end
        Hs{i,j} = D;
    end
end


% Kernel computation
fprintf(1,'\n');
disp('Computing kernel matrix...');
K = zeros(N,N);
fprintf(1,'%s','Progress: 0%');
for i=1:N
    
    last_percent = floor(100*(i-1)/N);
    cur_percent = floor(100*i/N);
    
    if last_percent ~= cur_percent
        last_str = sprintf('%d%s',last_percent,'%');
        fprintf(1,repmat('\b', 1, length(last_str)));
        fprintf(1,'%d%s',cur_percent,'%');
    end
    
    for j=i:N
        k = 0;
        intersec = zeros(L, 1);
        for p=1:L
            intersec(p) = sum(sum(min(Hs{i,p}, Hs{j,p})));
        end
        
        k = k + intersec(L);
        for p=1:L-1
            k = k + (1.0/(2^(L-p)))*(intersec(p)-intersec(p+1));
        end

        K(i,j) = k;
        K(j,i) = K(i,j);
    end
end
%runtime = toc;
fprintf(1,'\n');
%disp(['Kernel computation took ', num2str(runtime), ' sec']);
end