classdef utils_learned
    methods(Static)
        
        function theta_new = sample_hybrid_prior(N, naive_parts, naive_w, learned_lims)
            % Samples N particles.
            % Cols 1-4: Resampled from naive_parts (Bootstrap)
            % Cols 5-18: Uniform random within learned_lims, or exponential
            % distribution for the variability parameters (sigmaS's)
            
            D_learned = size(learned_lims, 1);
            
            % sample naive parameters
            % This preserves the correlation structure of the naive fit.
            idx = randsample(size(naive_parts,1), N, true, naive_w);
            theta_naive = naive_parts(idx, :);
            
            MU_C = 1;

            % sample learned parameters
            theta_learned = zeros(N, D_learned);
            for i = 1:D_learned % these are parameters 5-18, i.e indexed at 4 below the actual index
                low = learned_lims(i,1);
                high = learned_lims(i,2);
                if ismember(i,11:13) % 11, 12, 13 are sigmaS's 
                    % Truncated Exponential Distribution
                    pd = makedist('Exponential', 'mu', MU_C);
                    t_pd = truncate(pd, low, high);
                    theta_learned(:, i) = random(t_pd, N, 1);
                else
                    % Uniform Distribution
                    theta_learned(:,i) = low + (high - low) .* rand(N,1);
                end
            end
            
            theta_new = [theta_naive, theta_learned];
        end

        function lp = log_prior_hybrid(theta, limits, naive_mu, naive_Sigma)
            % Calculates Log Prior.
            % For Naive (1-4): Multivariate Normal approximation of previous posterior
            % For Learned (5-end): Uniform/exponential density.
            %
            % Inputs:
            %   theta: N x D matrix of parameters
            %   limits: D_learned x 2 matrix of uniform bounds
            %   naive_mu: 1 x 4 vector of means for naive params (Optional)
            %   naive_Sigma: 4 x 4 covariance matrix for naive params (Optional)

            [N, D_total] = size(theta);
            D_learned = size(limits, 1);
            D_naive = D_total - D_learned;

            % 1. Naive Part: Multivariate Normal Log PDF
            % We penalize drift away from the naive posterior region.
            if nargin >= 4 && ~isempty(naive_mu) && ~isempty(naive_Sigma)
                % Using mvnpdf for log probability density
                % theta(:, 1:4) are the naive parameters
                try
                    % mvnpdf returns the density, we take log.
                    % (mvnpdf supports vectorized input N x D)
                    p_naive = mvnpdf(theta(:, 1:D_naive), naive_mu, naive_Sigma);
                    lp_naive = log(p_naive);
                catch
                    % Fallback if Sigma is singular or issues arise
                    warning('Error calculating MVN prior. Defaulting to flat prior for Naive params.');
                    lp_naive = zeros(N, 1);
                end
            else
                % If moments aren't provided, we effectively use a flat prior
                % (Not recommended if you want to constrain the naive params)
                lp_naive = zeros(N, 1);
            end

            % 2. Learned Part: Uniform Log PDF
            lp_learned = zeros(N,1);

            for i = 1:D_learned
                col = D_naive + i;
                val = theta(:,col);
                low = limits(i,1);
                high = limits(i,2);

                in_bounds = (val >= low) & (val <= high);

                % If out of bounds, log prob is -Inf
                lp_learned(~in_bounds) = -Inf;


                if ismember(i,11:13) %i in 11, 12, 13 are sigmaS's
                    % Log PDF of Truncated Exponential
                    % f(x) = (lambda * exp(-lambda*x)) / (CDF(high) - CDF(low))
                    % log(f(x)) = log(lambda) - lambda*x - log(Normalization)

                    MU_C = 1; 
                    lambda = 1/MU_C;

                    % % Normalization constant Z = exp(-lam*low) - exp(-lam*high)
                    % % We factor out exp(-lam*low) for numerical stability if needed,
                    % % but standard calculation is usually fine for these ranges.
                    % cdf_diff = exp(-lambda * low) - exp(-lambda * high);
                    % log_Z = log(cdf_diff / lambda); % Integral of exp PDF over [low, high]

                    % Note: integral of lambda*exp(-lambda*x) is -exp(-lambda*x)
                    % Probability mass in range = [-exp(-lam*high) - (-exp(-lam*low))]
                    %                           = exp(-lam*low) - exp(-lam*high)

                    log_prob_mass = log(exp(-lambda*low) - exp(-lambda*high));

                    % Update log prior for valid particles
                    % log(lambda) - lambda*x - log(Mass)
                    lp_learned(in_bounds) = lp_learned(in_bounds) ...
                        + log(lambda) ...
                        - lambda * val(in_bounds) ...
                        - log_prob_mass;

                else % Uniform
                    % log(1 / (high - low))
                    lp_learned(in_bounds) = lp_learned(in_bounds) - log(high - low);
                end

                % % Add Log-Density of Uniform: log(1/(high-low)) = -log(range)
                % range = high - low;
                % lp_learned(in_bounds) = lp_learned(in_bounds) - log(range);
            end


            % Total Log Prior
            lp = lp_naive + lp_learned;

            % %%%%%% for mixture model only
            % p0pNegWrong = sum(theta(:,5:7) + theta(:,15:17) > 1, 2) > 0;
            % lp(p0pNegWrong) = -Inf;
        
        end

        function theta_new = sample_block_kernel(theta_old, cholC)
            % Hybrid Kernel with Rejuvenation.
            % Generates a standard Multivariate Normal perturbation based on cholC (the Cholesky Factor of the covariance matrix).
            % The code allows you to set the probability of mutating the
            % naive parameters on the idea that they shouldn't be updated
            % too much away from the naive posterior. However, in practice
            % we found that the algorithm gets stuck at singular values for
            % naive parameters if you don't mutate naive parameters more
            % often. Left the code like this to preserve the option of
            % reducing the probability of mutating naive parameters.
            
            REJUVENATION_PROB = 1; % chance to move naive params
            NAIVE_INDICES = 1:4;
            
            [N, D] = size(theta_old);
            
            % 1. Standard MVN Proposal
            perturbation = randn(N,D) * cholC; 
            
            % 2. Create Mask
            % mask = 1 if we Rejuvenate (Move All), 0 if we Block Naive (Move Learned Only)
            do_rejuv = rand(N,1) < REJUVENATION_PROB;
            
            % 3. Apply Mask
            % If NOT rejuvenating, set cols 1:4 of perturbation to 0
            perturbation(~do_rejuv, NAIVE_INDICES) = 0;
            
            theta_new = theta_old + perturbation;
        end


        function distResult = CvMdivergenceLearned(sim, obs)
            % sim and obs should be 7x4 cell arrays
            numDilutions = size(sim,1);
            numVoltages = size(sim,2);

            dists = zeros(numDilutions, numVoltages, 3, 3); % dilutions, voltages, metrics, epochs
            simCorrs = zeros(numDilutions, numVoltages, 3); % dilutions, voltages, metrics
            obsCorrs = zeros(numDilutions, numVoltages, 3); % dilutions, voltages, metrics
            for d=1:numDilutions
                for v=1:numVoltages
                    for i=1:3
                        for e=1:2
                            try
                                [~,~,dists(d,v,i,e)] = cmtest2(sim{d,v}(:,i,e), obs{d,v}(:,i,e));
                            catch
                                dists(d,v,i,e) = 1000; % if there are not enough non-Nan samples to calculate CvM distance
                            end
                        end
                        try
                            [~,~,dists(d,v,i,3)] = cmtest2(sim{d,v}(:,i,2)-sim{d,v}(:,i,1), obs{d,v}(:,i,2)-obs{d,v}(:,i,1));
                        catch
                            dists(d,v,i,3) = 1000; % if there are not enough non-Nan samples to calculate CvM distance
                        end
                        try
                            nans = isnan(sim{d,v}(:,i,1)) | isnan(sim{d,v}(:,i,2));
                            simCorrs(d,v,i) = corr(sim{d,v}(~nans,i,1), sim{d,v}(~nans,i,2),'Type','Kendall');
                            nans = isnan(obs{d,v}(:,i,1)) | isnan(obs{d,v}(:,i,2));
                            obsCorrs(d,v,i) = corr(obs{d,v}(~nans,i,1), obs{d,v}(~nans,i,2),'Type','Kendall');
                            
                        catch
                        end
                    end
                end
            end
            corrDiffs = squeeze(mean(simCorrs - obsCorrs,[1 2]));
            distResult = norm(dists(:)) + 6*norm(corrDiffs); % epsilon is around 5 for CvM dists; corr diffs around 0.17. Times 6 would make the corr diffs weighted around 0.2.
        end


    end
end