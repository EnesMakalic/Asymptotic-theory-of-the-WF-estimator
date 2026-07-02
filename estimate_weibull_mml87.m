function [k_mml, lambda_mml, codelength, exitflag] = estimate_weibull_mml87(x, options)
% ESTIMATE_WEIBULL_MML87 Efficient MML87 estimation for Weibull distribution
%
% [k_mml, lambda_mml, codelength, exitflag] = estimate_weibull_mml87(x)
% [k_mml, lambda_mml, codelength, exitflag] = estimate_weibull_mml87(x, Name=Value)
%
% Inputs:
%   x - vector of positive observations (n x 1 or 1 x n)
%
% Optional Name-Value pairs:
%   InitMethod - 'moments' (default), 'mle', or 'custom'
%   InitParams - [k0, lambda0] if InitMethod is 'custom'
%   Display    - 'off' (default), 'iter', 'final'
%
% Outputs:
%   k_mml      - MML87 estimate of shape parameter
%   lambda_mml - MML87 estimate of scale parameter  
%   codelength - minimum MML87 codelength (in nats)
%   exitflag   - optimization exit status (1 = success)
%
% Example:
%   x = wblrnd(1.5, 2.5, 100, 1);
%   [k, lambda] = estimate_weibull_mml87(x);
%   [k, lambda] = estimate_weibull_mml87(x, InitMethod='mle', Display='iter');
arguments
    x (:,1) double {mustBePositive, mustBeNonempty};
    options.InitMethod (1,1) string {mustBeMember(options.InitMethod, ...
        ["moments", "mle", "custom"])} = "moments";
    options.InitParams (1,2) double {mustBePositive} = [1,1];
    options.Display (1,1) string {mustBeMember(options.Display, ...
        ["off", "iter", "final"])} = "off";
end

n = length(x);

% Pre-compute constants (persistent for efficiency)
persistent KAPPA2
if isempty(KAPPA2)
    KAPPA2 = 5 / (36 * sqrt(3));
end

% Pre-compute data statistics for efficiency
log_x = log(x);
sum_log_x = sum(log_x);

% Get initial parameters
theta0 = get_initial_params(x, options.InitMethod, options.InitParams);

% Set up optimization objective
objective = @(log_theta) mml87_objective(log_theta, x, log_x, sum_log_x, n, KAPPA2);

% Optimization options
opt_options = optimoptions('fminunc', ...
    'Display', char(options.Display), ...
    'Algorithm', 'quasi-newton', ...
    'SpecifyObjectiveGradient', true, ...
    'MaxFunctionEvaluations', 1000, ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-10);

% Optimize in log-space for unconstrained optimization
log_theta0 = log(theta0);
[log_theta_opt, codelength, exitflag] = fminunc(objective, log_theta0, opt_options);

% Transform back to original parameters
k_mml = exp(log_theta_opt(1));
lambda_mml = exp(log_theta_opt(2));

% Handle convergence warning
if exitflag <= 0
    warning('Optimization did not converge properly. Results may be suboptimal.');
end
end

function theta0 = get_initial_params(x, method, custom_params)
% GET_INITIAL_PARAMS Compute initial parameter estimates

arguments
    x (:,1) double
    method (1,1) string
    custom_params (1,2) double
end

switch method
    case "custom"
        if isempty(custom_params)
            error('InitParams must be provided when InitMethod is ''custom''');
        end
        theta0 = custom_params;
        
    case "mle"
        % Use MATLAB's built-in MLE (if Statistics Toolbox available)
        try
            phat = wblfit(x);
            theta0 = [phat(2), phat(1)];  % Note: wblfit returns [lambda, k]
        catch
            % Fall back to moments if wblfit not available
            theta0 = get_initial_params(x, "moments", []);
        end
        
    case "moments"
        % Method of moments initialization
        x_mean = mean(x);
        x_std = std(x);
        cv = x_std / x_mean;  % Coefficient of variation
        
        % Solve for k using CV relationship
        % CV = sqrt(Γ(1+2/k)/Γ(1+1/k)^2 - 1)
        % Use approximation for initial guess
        if cv < 1
            k0 = 1.086 / cv;  % Empirical approximation
        else
            k0 = 1.0;
        end
        
        % Ensure reasonable range
        k0 = max(0.5, min(k0, 10));
        
        % Compute lambda from mean
        lambda0 = x_mean / gamma(1 + 1/k0);
        
        theta0 = [k0, lambda0];
end
end

function [f, g] = mml87_objective(log_theta, x, log_x, sum_log_x, n, KAPPA2)
% MML87_OBJECTIVE Compute objective and gradient efficiently
%
% Uses log-transformed parameters for unconstrained optimization
% Fisher determinant: det(J) = n^2 * pi^2 / (6 * lambda^2)

arguments
    log_theta (2,1) double
    x (:,1) double
    log_x (:,1) double
    sum_log_x (1,1) double
    n (1,1) double
    KAPPA2 (1,1) double
end

% Transform parameters back
k = exp(log_theta(1));
lambda = exp(log_theta(2));

% Precompute common terms
log_lambda = log_theta(2);
x_over_lambda = x / lambda;
x_over_lambda_k = x_over_lambda.^k;
log_x_over_lambda = log_x - log_lambda;

% Log-likelihood
log_lik = n*log(k) - n*log_lambda + (k-1)*sum(log_x_over_lambda) - sum(x_over_lambda_k);

% Half-Cauchy prior
log_prior_k = log(2) - log(pi) - log(1 + k^2);
log_prior_lambda = log(2) - log(pi) - log(1 + lambda^2);
log_prior = log_prior_k + log_prior_lambda;

% Log determinant of Fisher information matrix
% det(J) = n^2 * pi^2 / (6 * lambda^2)
log_det_J = 2*log(n) + 2*log(pi) - log(6) - 2*log_lambda;

% MML87 codelength
assertion = -log_prior + 0.5*log_det_J + log(KAPPA2);
detail = -log_lik + 1;  % p/2 = 2/2 = 1
f = assertion + detail;

% Compute gradient if requested
if nargout > 1
    % Gradient w.r.t. k
    dloglik_dk = n/k + sum(log_x_over_lambda) - sum(x_over_lambda_k .* log_x_over_lambda);
    dlogprior_dk = -2*k / (1 + k^2);
    
    % Gradient of log(det(J)) w.r.t. k
    % Since det(J) doesn't depend on k, derivative is 0
    dlogdet_dk = 0;
    
    df_dk = -dlogprior_dk + 0.5*dlogdet_dk - dloglik_dk;
    
    % Gradient w.r.t. lambda
    dloglik_dlambda = -n/lambda - (k-1)*n/lambda + k*sum(x_over_lambda_k)/lambda;
    dlogprior_dlambda = -2*lambda / (1 + lambda^2);
    
    % Gradient of log(det(J)) w.r.t. lambda
    % d/dlambda[2*log(n) + 2*log(pi) - log(6) - 2*log(lambda)] = -2/lambda
    dlogdet_dlambda = -2/lambda;
    
    df_dlambda = -dlogprior_dlambda + 0.5*dlogdet_dlambda - dloglik_dlambda;
    
    % Chain rule for log-transformed parameters
    g = [df_dk * k; df_dlambda * lambda];
end

end
