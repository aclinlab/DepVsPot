classdef utils_naive_combined_MS
    methods(Static)
        %% prior sampler

        function theta = priorSamplerNaiveFn(N, mu, b)
            % Draw N parameter vectors from prior
            % abc_smc_adaptive expects theta to be N x D
            % ie each particle is one row; columns are parameters

            D = 5;
            % theta = zeros(N, D);
            % theta will be in this order:
            % 1. muXenter
            % 2. sigmaXenter
            % 3. slope_muXstay
            % 4. slope_sigmaXstay
            % 5. model number

            limits = ...
                [-2 2;
                0 2;
                -5 10;
                -5 10];
            limitRange = diff(limits,1,2);
            if size(limits,1) ~= D-1
                error('number of dimensions in prior limits does not match')
            end
            theta = nan(N,D); % inactive parameters default to NaN as error-catching measure

            theta(:,5) = randi(3,N,1); % assign model numbers

            theta(:,1) = rand(N,1).*repmat(limitRange(1)',N,1) + repmat(limits(1,1)',N,1);

            % sample from a Laplace distribution centered on 0
            u = rand(N,1) - 0.5;
            theta(:,3) = -b.*sign(u) .* log(1-2*abs(u));

            pd = makedist('Exponential','mu',mu); 
            % if model 2 or 3, fill theta(2) with prior
            idx23 = ismember(theta(:,5), [2 3]); % true(N,1); % use this for negative control with identical models
            truncated_pd = truncate(pd,limits(2,1),limits(2,2));
            theta(idx23,2) = random(truncated_pd,sum(idx23),1); % get N=sum(idx23) particles

            % if model 3, fill theta(4)
            idx3 = theta(:,5) == 3; % true(N,1); % use this for negative control with identical models
            numModel3 = sum(idx3);
            u = rand(numModel3,1) - 0.5;
            theta(idx3,4) = -b.*sign(u) .* log(1-2*abs(u));

        end


        %% prior log prob

        function lp = priorLogPdfNaiveFn(theta, mu, b)
            % Log prior pdf for NxD theta
            % Must return -Inf for impossible values
            % constants are removed because this is only used for prior ratio
            % calculations - MUST ADD THEM BACK IN IF USED FOR MODEL SELECTION!
            [N,~] = size(theta);
            lp = zeros(N,1);
            limits = ...
                [-2 2;
                0 2;
                -5 10;
                -5 10];

            in = true(N,1);

            % check limits for active params only
            for i=[1 3] % uncomments for different models % 1:4 use this for negative control with identical models
                in = in & theta(:,i) >= limits(i,1) & theta(:,i) <= limits(i,2);
            end

            % model index (m) must be valid
            m = round(theta(:,5));
            in = in & ismember(m, [1 2 3]);

            lambda = 1/mu;

            % uncomment this part for different models

            % Laplace prior for slope_mu_stay
            lp = lp - abs(theta(:,3))/b;

            % sigmaXenter active in models 2 and 3
            idx23 = in & (m == 2 | m == 3);
            lp(idx23) = lp(idx23) - lambda * theta(idx23, 2);

            % slope_sigma_stay active in model 3 only. Laplace prior for slope_sigma_stay
            idx3 = in & (m == 3);
            lp(idx3) = lp(idx3) - abs(theta(idx3, 4))/b;

            % bounds check on sigma params for active models
            idx23_out_of_bounds = idx23 & ~(theta(:,2) >= limits(2,1) & theta(:,2) <= limits(2,2));
            in(idx23_out_of_bounds) = false;
            idx3_out_of_bounds = idx3 & ~(theta(:,4) >= limits(4,1) & theta(:,4) <= limits(4,2));
            in(idx3_out_of_bounds) = false;

            % 
            % 
            % % use this part for negative control for identical models
            % % this is the equation for log prior of exponential distribution
            % lp = lp  - lambda*(theta(:,2)) ... % - log(diff(limits(1,:))) + log(lambda) % constants
            %     - abs(theta(:,3))/b ... % log prior of Laplace
            %     - abs(theta(:,4))/b; % log prior of Laplace% - log(2*pi*sigma^2) % constants
            %     % - 1/(2*sigma^2)*theta(:,3) ... % and this is the log prior of Gaussian 
            %     % - 1/(2*sigma^2)*(theta(:,3).^2 + theta(:,4).^2); % and this is the log prior of Gaussian % - log(2*pi*sigma^2) % constants

            lp(~in) = -Inf;

        end



        %% distance

        % To test: whether it helps to base the distance only on P(enter)
        % and P(stay) (not P(approach))

        function cmdist = CvMdivergenceNaive(sim, obs)
            % sim and obs should be 7x1 cell arrays
            numDilutions = length(sim);
            dists = zeros(numDilutions, 3);
            for d=1:numDilutions
                for i=1:3
                    [~,~,dists(d,i)] = cmtest2(sim{d}(:,i), obs{d}(:,i));
                end
            end
            cmdist = norm(dists(:));
        end


        function theta_star = kernelSamplerFn(theta_old, cholC)
            % cholC is either:
            %   - a matrix (standard case, no model selection)
            %   - a cell array of per-model matrices (model selection case)
            % Symmetric Gaussian RW kernel
            % theta_old: 1×D
            % cholC: D×D cholesky decomposition of covariance matrix
            % Returns 1×D vector

            theta_star = theta_old;  % copy, preserving NaNs in inactive cols

            activeCols = {[1 3], [1 2 3], [1 2 3 4]}; % {1:4, 1:4, 1:4}; % for negative control comparing identical models
            if iscell(cholC)
                % model selection case
                m = round(theta_old(end));  % model index is last column
                cholC_m = cholC{m};
                ac = activeCols{m};  % get active columns for this model
                D_m = length(ac);
                theta_star(ac) = theta_old(ac) + randn(1, D_m) * cholC_m;
            else
                % standard case: perturb all columns
                D = length(theta_old);
                theta_star = theta_old + randn(1, D) * cholC; % Del Moral et al. suggest 2*C but we think 1*C is better?
            end
        end



    end
end
