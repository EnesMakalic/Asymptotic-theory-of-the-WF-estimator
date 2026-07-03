function [fig1, fig2] = plot_weibull_bias(niter, options)
% PLOT_WEIBULL_BIAS Plot empirical bias of MML and MLE estimators for Weibull k
%
% [fig1, fig2] = plot_weibull_bias(niter)
% [fig1, fig2] = plot_weibull_bias(niter, Name=Value)
%
% Inputs:
%   niter - number of Monte Carlo iterations for bias estimation
%
% Optional Name-Value pairs:
%   k_range     - [k_min, k_max] for range of k values (default: [0.1, 10])
%   n_k_values  - number of k values to test (default: 20)
%   sample_size - sample size for each iteration (default: 20)
%   lambda      - fixed scale parameter (default: 1)
%   show_stderr - logical, show standard error bands (default: false)
%
% Outputs:
%   fig1 - figure handle for bias comparison plot
%   fig2 - figure handle for empirical vs analytical bias plot
%
% Example:
%   plot_weibull_bias(1000);
%   plot_weibull_bias(500, k_range=[0.5, 5], sample_size=50);
%   [f1, f2] = plot_weibull_bias(1000, show_stderr=true);

arguments
    niter (1,1) double {mustBePositive, mustBeInteger}
    options.k_range (1,2) double {mustBePositive} = [0.1, 10]
    options.n_k_values (1,1) double {mustBePositive, mustBeInteger} = 20
    options.sample_size (1,1) double {mustBePositive, mustBeInteger} = 20
    options.lambda (1,1) double {mustBePositive} = 1
    options.show_stderr (1,1) logical = false
end

% Validate k_range
if options.k_range(1) >= options.k_range(2)
    error('k_range(1) must be less than k_range(2)');
end

% Generate k values (log-spaced for better coverage)
k_true = logspace(log10(options.k_range(1)), log10(options.k_range(2)), ...
                  options.n_k_values);

% Pre-allocate arrays
bias_mml = zeros(size(k_true));
bias_mle = zeros(size(k_true));
stderr_mml = zeros(size(k_true));
stderr_mle = zeros(size(k_true));

% Analytical bias formulas
bias_mle_analytical = zeros(size(k_true));
bias_mml_analytical = zeros(size(k_true));

% Progress display
fprintf('Computing empirical bias for %d k values with %d iterations each...\n', ...
        options.n_k_values, niter);
fprintf('Sample size: n = %d, lambda = %.2f\n\n', options.sample_size, options.lambda);

% Start timer
tic;

% Loop over true k values
for i = 1:length(k_true)
    k_val = k_true(i);
    
    % Compute analytical bias
    [bias_mle_analytical(i), bias_mml_analytical(i)] = ...
        analytical_bias_weibull(k_val, options.lambda, options.sample_size);
    
    % Storage for estimates
    k_mml_estimates = zeros(niter, 1);
    k_mle_estimates = zeros(niter, 1);
    
    % Monte Carlo simulation
    parfor iter = 1:niter
        % Generate data
        x = wblrnd(options.lambda, k_val, options.sample_size, 1);
        
        % MML estimate (suppress warnings)
        try
            [k_mml_estimates(iter), ~] = estimate_weibull_mml87(x, Display="off",InitMethod="mle");
        catch
            k_mml_estimates(iter) = NaN;
        end
        
        % MLE estimate
        try
            phat = wblfit(x);
            k_mle_estimates(iter) = phat(2);
        catch
            k_mle_estimates(iter) = NaN;
        end
    end
    
    % Remove NaN values
    k_mml_valid = k_mml_estimates(~isnan(k_mml_estimates));
    k_mle_valid = k_mle_estimates(~isnan(k_mle_estimates));
    
    % Compute empirical bias
    bias_mml(i) = mean(k_mml_valid) - k_val;
    bias_mle(i) = mean(k_mle_valid) - k_val;
    
    % Compute standard error
    stderr_mml(i) = std(k_mml_valid) / sqrt(length(k_mml_valid));
    stderr_mle(i) = std(k_mle_valid) / sqrt(length(k_mle_valid));
    
    % Progress indicator
    if mod(i, max(1, floor(options.n_k_values/10))) == 0
        fprintf('Progress: %d/%d (%.0f%%) - Elapsed: %.1f sec\n', ...
                i, options.n_k_values, 100*i/options.n_k_values, toc);
    end
end

elapsed_time = toc;
fprintf('Done! Total time: %.1f seconds\n\n', elapsed_time);

% =========================================================================
% FIGURE 1: Empirical MML vs MLE Bias
% =========================================================================
fig1 = figure('Position', [100, 100, 900, 600]);
hold on; grid on;

% Plot standard error bands if requested
if options.show_stderr
    % MML error band
    fill([k_true, fliplr(k_true)], ...
         [bias_mml + 2*stderr_mml, fliplr(bias_mml - 2*stderr_mml)], ...
         'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'DisplayName', 'WF ±2 SE');
    
    % MLE error band
    fill([k_true, fliplr(k_true)], ...
         [bias_mle + 2*stderr_mle, fliplr(bias_mle - 2*stderr_mle)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
         'DisplayName', 'MLE ±2 SE');
end

% Plot bias curves
% plot(k_true, bias_mml, 'b-', 'LineWidth', 2.5, ...
%      'DisplayName', 'WF Bias', 'Marker', 'o', 'MarkerSize', 6);
plot(k_true, bias_mml, 'k-', 'LineWidth', 2.0, ...
     'DisplayName', 'WF Bias', 'Marker', 'o', 'MarkerSize', 6);
% plot(k_true, bias_mle, 'r--', 'LineWidth', 2.0, ...
%      'DisplayName', 'MLE Bias', 'Marker', 's', 'MarkerSize', 6);
plot(k_true, bias_mle, 'k--', 'LineWidth', 2.0, ...
     'DisplayName', 'MLE Bias', 'Marker', 's', 'MarkerSize', 6);
%plot(k_true, zeros(size(k_true)), 'k:', 'LineWidth', 1.5, ...
%     'HandleVisibility', 'off');

% Labels and formatting
xlabel('True Shape Parameter (k)', 'FontSize', 18);
ylabel('Empirical Bias', 'FontSize', 18);
%title(sprintf('Weibull Shape Parameter Bias: WF vs MLE\nn=%d, \\lambda=%.0f', ...
%              options.sample_size, options.lambda), ...
%      'FontSize', 18);
legend('Location', 'northwest', 'FontSize', 18);

% Set x-axis to log scale if range spans more than one order of magnitude
if options.k_range(2) / options.k_range(1) > 10
    set(gca, 'XScale', 'log');
    xlabel('Shape Parameter (k)', 'FontSize', 18);
end

% Improve appearance
set(gca, 'FontSize', 18, 'LineWidth', 1);
set(gcf,'Color','w');
box on;
%yline(0, 'k-', 'LineWidth', 1);

% =========================================================================
% FIGURE 2: Empirical vs Analytical Bias
% =========================================================================
fig2 = figure('Position', [150, 150, 1400, 600]);

% Left subplot: MLE bias comparison
subplot(1, 2, 1);
hold on; grid on;

% Plot empirical vs analytical for MLE
%plot(k_true, bias_mle, 'ro-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
%     'DisplayName', 'Empirical MLE Bias', 'MarkerFaceColor', 'r');
plot(k_true, bias_mle, 'ko-', 'LineWidth', 2.0, 'MarkerSize', 8, ...
     'DisplayName', 'Empirical MLE Bias', 'MarkerFaceColor', 'k');
%plot(k_true, bias_mle_analytical, 'b--', 'LineWidth', 2.5, ...
%     'DisplayName', 'Analytical MLE Bias (O(n^{-1}))', 'Marker', 's', 'MarkerSize', 6);
plot(k_true, bias_mle_analytical, 'k--', 'LineWidth', 2.0, ...
     'DisplayName', 'Analytical MLE Bias (O(n^{-1}))', 'Marker', 's', 'MarkerSize', 8);
%yline(0, 'k:', 'LineWidth', 1.5);

xlabel('Shape Parameter (k)', 'FontSize', 18);
ylabel('Bias', 'FontSize', 18);
%title(sprintf('MLE Bias: Empirical vs Analytical\nn=%d, \\lambda=%.0f', ...
%              options.sample_size, options.lambda), ...
%      'FontSize', 18);
legend('Location', 'northwest', 'FontSize', 18);

if options.k_range(2) / options.k_range(1) > 10
    set(gca, 'XScale', 'log');
end
set(gca, 'FontSize', 18, 'LineWidth', 1);
box on;

% Right subplot: MML bias comparison
subplot(1, 2, 2);
hold on; grid on;

% Plot empirical vs analytical for MML
% plot(k_true, bias_mml, 'go-', 'LineWidth', 2.5, 'MarkerSize', 8, ...
%      'DisplayName', 'Empirical MML Bias', 'MarkerFaceColor', 'g');
plot(k_true, bias_mml, 'ko-', 'LineWidth', 2.0, 'MarkerSize', 8, ...
     'DisplayName', 'Empirical WF Bias', 'MarkerFaceColor', 'k');
% plot(k_true, bias_mml_analytical, 'm--', 'LineWidth', 2.5, ...
%      'DisplayName', 'Analytical MML Bias (O(n^{-1}))', 'Marker', 'd', 'MarkerSize', 6);
plot(k_true, bias_mml_analytical, 'k--', 'LineWidth', 2.0, ...
     'DisplayName', 'Analytical WF Bias (O(n^{-1}))', 'Marker', 'd', 'MarkerSize', 8);
%yline(0, 'k:', 'LineWidth', 1.5);

xlabel('Shape Parameter (k)', 'FontSize', 18);
ylabel('Bias', 'FontSize', 18);
%title(sprintf('WF Bias: Empirical vs Analytical\nn=%d, \\lambda=%.0f', ...
%              options.sample_size, options.lambda), ...
%      'FontSize', 18);
legend('Location', 'northwest', 'FontSize', 18);

if options.k_range(2) / options.k_range(1) > 10
    set(gca, 'XScale', 'log');
end
set(gca, 'FontSize', 18, 'LineWidth', 1);
set(gcf,'Color','w');
box on;

% =========================================================================
% Print summary statistics
% =========================================================================
fprintf('Summary Statistics:\n');
fprintf('─────────────────────────────────────────────────────────────────\n');
fprintf('  MLE:\n');
fprintf('    Mean absolute empirical bias:  %.6f\n', mean(abs(bias_mle)));
fprintf('    Mean absolute analytical bias: %.6f\n', mean(abs(bias_mle_analytical)));
fprintf('    Max absolute empirical bias:   %.6f (at k=%.2f)\n', ...
        max(abs(bias_mle)), k_true(abs(bias_mle)==max(abs(bias_mle))));
fprintf('    RMS empirical bias:            %.6f\n', sqrt(mean(bias_mle.^2)));
fprintf('    RMS analytical bias:           %.6f\n', sqrt(mean(bias_mle_analytical.^2)));

% Compute correlation between empirical and analytical
corr_mle = corr(bias_mle(:), bias_mle_analytical(:));
fprintf('    Correlation (empirical-analytical): %.4f\n', corr_mle);
fprintf('\n');

fprintf('  WF:\n');
fprintf('    Mean absolute empirical bias:  %.6f\n', mean(abs(bias_mml)));
fprintf('    Mean absolute analytical bias: %.6f\n', mean(abs(bias_mml_analytical)));
fprintf('    Max absolute empirical bias:   %.6f (at k=%.2f)\n', ...
        max(abs(bias_mml)), k_true(abs(bias_mml)==max(abs(bias_mml))));
fprintf('    RMS empirical bias:            %.6f\n', sqrt(mean(bias_mml.^2)));
fprintf('    RMS analytical bias:           %.6f\n', sqrt(mean(bias_mml_analytical.^2)));

% Compute correlation between empirical and analytical
corr_mml = corr(bias_mml(:), bias_mml_analytical(:));
fprintf('    Correlation (empirical-analytical): %.4f\n', corr_mml);
fprintf('\n');

fprintf('  Bias Reduction:\n');
if mean(abs(bias_mml)) < mean(abs(bias_mle))
    fprintf('    WF has %.1f%% lower average absolute bias than MLE\n', ...
            100*(1 - mean(abs(bias_mml))/mean(abs(bias_mle))));
else
    fprintf('    MLE has %.1f%% lower average absolute bias than WF\n', ...
            100*(1 - mean(abs(bias_mle))/mean(abs(bias_mml))));
end

% Compare analytical predictions
fprintf('    Analytical: WF has %.1f%% lower bias than MLE\n', ...
        100*(1 - mean(abs(bias_mml_analytical))/mean(abs(bias_mle_analytical))));
fprintf('─────────────────────────────────────────────────────────────────\n');

% Save results to base workspace
results = struct('k_true', k_true, ...
                 'bias_mml', bias_mml, 'bias_mle', bias_mle, ...
                 'bias_mml_analytical', bias_mml_analytical, ...
                 'bias_mle_analytical', bias_mle_analytical, ...
                 'stderr_mml', stderr_mml, 'stderr_mle', stderr_mle, ...
                 'niter', niter, 'sample_size', options.sample_size, ...
                 'lambda', options.lambda);
assignin('base', 'bias_results', results);
fprintf('\nResults saved to workspace variable ''bias_results''\n');

end

function [bias_mle_k, bias_mml_k] = analytical_bias_weibull(k, lambda, n)
% ANALYTICAL_BIAS_WEIBULL Compute analytical bias formulas for Weibull distribution
%
%
% Inputs:
%   k      - true shape parameter
%   lambda - true scale parameter
%   n      - sample size
%
% Outputs:
%   bias_mle_k - analytical MLE bias for shape parameter k
%   bias_mml_k - analytical MML bias for shape parameter k

arguments
    k (1,1) double {mustBePositive}
    lambda (1,1) double {mustBePositive}
    n (1,1) double {mustBePositive}
end

% Mathematical constants
gamma_const = 0.577215664901532860606512090082;  % Euler-Mascheroni constant
zeta3 = 1.202056903159594285399738161511;        % Riemann zeta(3)

% MLE bias for k 
bias_mle_k = (18 * k * (pi^2 - 2*zeta3)) / (n * pi^4);

% MML bias for k 
% First component of the shift vector
correction_k = (6 / pi^2) * ((gamma_const - 1) * (lambda^2 - 1) / (lambda^2 + 1) ...
                             - (2 * k^3) / (k^2 + 1));

bias_mml_k = bias_mle_k + correction_k / n;
end
