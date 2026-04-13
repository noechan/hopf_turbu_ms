function precompute_hopf_geometry(data_path, output_path, lambda)

    % Load coordinates
    cog_file = fullfile(data_path, 'SchaeferCOG.mat');
    if ~isfile(cog_file)
        error('SchaeferCOG.mat not found at: %s', cog_file);
    end
    load(cog_file, 'SchaeferCOG');

    NPARCELLS = size(SchaeferCOG,1);
    NR = 400;

    % Optional: define output file
    outfile = fullfile(output_path, 'precomputed_geometry.mat');

    % Create output folder if needed
    if ~isfolder(output_path)
        mkdir(output_path);
    end

    % Compute rr exactly as in original scripts
    rr = zeros(NPARCELLS, NPARCELLS);
    for i = 1:NPARCELLS
        for j = 1:NPARCELLS
            rr(i,j) = norm(SchaeferCOG(i,:) - SchaeferCOG(j,:));
        end
    end

    range = max(max(rr));
    delta = range / NR;

    xrange = zeros(1, NR);
    for i = 1:NR
        xrange(i) = delta/2 + delta*(i-1);
    end

    C = zeros(NPARCELLS, NPARCELLS);

    LAMBDA = [0.27 0.24 0.21 0.18 0.15 0.12 0.09 0.06 0.03 0.01];

    NLAMBDA = length(LAMBDA);
    C1 = zeros(NLAMBDA, NPARCELLS, NPARCELLS);
    [~, indsca] = min(abs(LAMBDA - lambda));

    ilam = 1;
    for lambda2 = LAMBDA
        for i = 1:NPARCELLS
            for j = 1:NPARCELLS
                C1(ilam,i,j) = exp(-lambda2 * rr(i,j));
            end
        end
        ilam = ilam + 1;
    end

    save(outfile, 'rr', 'range', 'delta', 'xrange', 'C', 'LAMBDA', 'NLAMBDA', 'C1', 'indsca', '-v7.3');

    fprintf('Precomputed geometry saved to:\n%s\n', outfile);

end