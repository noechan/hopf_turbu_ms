function fit = fit_gec_from_targets(C_init, FCemp, COVtauemp, f_diff, cfg)
%FIT_GEC_FROM_TARGETS Fit group-level generative effective connectivity.
%
% Inputs
%   C_init    : initial SC matrix (already scaled)
%   FCemp     : empirical FC target
%   COVtauemp : empirical lagged covariance target (normalised)
%   f_diff    : frequency difference input for hopf_int
%   cfg       : configuration struct from get_gec_config()
%
% Output
%   fit       : struct with fitted matrices and optimisation history

N = cfg.N;
Tau = cfg.Tau;
TR = cfg.TR;

C = C_init;
Cnew = C_init;
olderror = 100000;

errorFC = nan(cfg.maxIter, 1);
errorCOVtau = nan(cfg.maxIter, 1);

for iter = 1:cfg.maxIter

    % Linear Hopf approximation
    [FCsim, COVsim, COVsimtotal, A] = hopf_int(Cnew, f_diff(:,:), cfg.sigma);

    % Lagged covariance prediction
    COVtausim = expm((Tau * TR) * A) * COVsimtotal;
    COVtausim = COVtausim(1:N, 1:N);

    % Variance normalisation, preserving original implementation
    sigratiosim = zeros(N, N);
    for i = 1:N
        for j = 1:N
            sigratiosim(i,j) = 1 / sqrt(COVsim(i,i)) / sqrt(COVsim(j,j));
        end
    end
    COVtausim = COVtausim .* sigratiosim;

    % Error terms
    errorFC(iter) = mean(mean((FCemp - FCsim).^2));
    errorCOVtau(iter) = mean(mean((COVtauemp - COVtausim).^2));

    % Stopping rule, preserving original logic
    if mod(iter, cfg.stop_check_every) < 0.1
        errornow = mean(mean((FCemp - FCsim).^2)) + ...
                   mean(mean((COVtauemp - COVtausim).^2));

        if (olderror - errornow) / errornow < cfg.stop_rel_improvement
            break;
        end

        if olderror < errornow
            break;
        end

        olderror = errornow;
    end

    % Gradient-descent-like learning
    for i = 1:N
        for j = 1:N
            if (C(i,j) > 0 || i == j + cfg.homologue_offset || j == i + cfg.homologue_offset)
                Cnew(i,j) = Cnew(i,j) ...
                    + cfg.epsFC    * (FCemp(i,j)    - FCsim(i,j)) ...
                    + cfg.epsFCtau * (COVtauemp(i,j) - COVtausim(i,j));

                if Cnew(i,j) < 0
                    Cnew(i,j) = 0;
                end
            end
        end
    end

    Cnew = Cnew / max(max(Cnew)) * cfg.maxC;
end

last_iter = iter;

fit = struct();
fit.Ceffgroup = Cnew;
fit.FCsim = FCsim;
fit.COVsim = COVsim;
fit.COVtausim = COVtausim;
fit.A = A;
fit.errorFC = errorFC(1:last_iter);
fit.errorCOVtau = errorCOVtau(1:last_iter);
fit.n_iter = last_iter;
fit.final_error = errorFC(last_iter) + errorCOVtau(last_iter);
end