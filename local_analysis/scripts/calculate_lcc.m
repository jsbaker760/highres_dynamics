function lcc_all = calculate_lcc(DM,stepsize)
% Compute threshold values.
EdgeLengths = 0:stepsize:max(DM(:));
S = size(DM,1);        % DM is SxS (symmetric)
E = numel(EdgeLengths);

% Preallocate results (single precision for speed)
LCC_all_isolates = zeros(E, S, 'single');
E_all_isolates = repmat(EdgeLengths', 1, S);

% Set up a lightweight progress indicator for the parfor loop.
p = 0;  % progress counter
dq = parallel.pool.DataQueue;
afterEach(dq, @updateProgress);

% Iterate over each threshold value in parallel.
parfor i = 1:E
    e = EdgeLengths(i);
    % Create binary adjacency matrix (in single precision)
    G = single(DM <= e);
    % Remove self-loops with fast linear indexing
    G(1:S+1:end) = 0;

    % Compute node degrees (sum along rows)
    d = sum(G, 2);

    % Compute diagonal of G^3:
    % Note: For an undirected graph, (G^3)_ii = 2*(# of triangles at node i)
    tri = sum((G * G) .* G, 2);

    % Calculate local clustering coefficient:
    LCC = zeros(S, 1, 'single');
    valid = d >= 2;  % Only nodes with degree >= 2 have a defined clustering coefficient
    LCC(valid) = tri(valid) ./ (d(valid) .* (d(valid) - 1));

    LCC_all_isolates(i, :) = LCC;
    % Send a progress update.
    send(dq, 1);
end

% Return results in a structure.
lcc_all = struct;
lcc_all.LCC_all_isolates = LCC_all_isolates;
lcc_all.E_all_isolates = E_all_isolates;
% Nested function to update progress.
    function updateProgress(~)
        p = p + 1;
        fprintf('Completed %d of %d iterations\n', p, E);
    end
end
