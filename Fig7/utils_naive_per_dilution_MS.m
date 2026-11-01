classdef utils_naive_per_dilution_MS
    methods(Static)


        function result = uniformPriorLimits()

            result = [-2 2;
                0 2;
                -2 2;
                0 2];

        end

        function theta = priorSamplerFn(N, mu)
            % Draw N parameter vectors from prior
            % abc_smc_adaptive expects theta to be N x D
            % ie each particle is one row; columns are parameters

            D = 5;
            % theta = zeros(N, D);
            % theta will be in this order:
            % 1. muXenter
            % 2. sigmaXenter
            % 3. muXstay
            % 4. sigmaXstay

            limits = uniformPriorLimits();
            limitRange = diff(limits,1,2);
            if size(limits,1) ~= 4 %D
                error('number of dimensions in prior limits does not match')
            end
            theta = nan(N,D); % inactive parameters default to NaN as error-catching measure
            theta(:,[1 3]) = rand(N,2).*repmat(limitRange([1 3])',N,1) + repmat(limits([1 3],1)',N,1);

            theta(:,5) = randi(2,N,1);

            pd = makedist('Exponential','mu',mu);
            
            % if model 2, fill theta(2) and theta(4) with prior
            idx2 = theta(:,5) == 2;
            for i=[2 4]
                truncated_pd = truncate(pd,limits(i,1),limits(i,2));
                theta(idx2,i) = random(truncated_pd,sum(idx2),1);
            end
        end

        %% prior log prob

        function lp = priorLogPdfFn(theta, mu)
            % Log prior pdf for NxD theta
            %
            % Must return -Inf for impossible values
            [N,~] = size(theta);
            lp = zeros(N,1);
            limits = uniformPriorLimits(); % 4x2, only continuous params

            % bounds check on active params only
            in = theta(:,1) >= limits(1,1) & theta(:,1) <= limits(1,2) & ...
                theta(:,3) >= limits(3,1) & theta(:,3) <= limits(3,2);

            % model index (m) must be valid
            m = round(theta(:,5));
            % in = in & ismember(m, [1 2 3]);
            in = in & ismember(m, [1 2]);

            % log priors for exponential parameters
            % mu = 0.05; %take as a parameter instead
            lambda = 1/mu;

            % sigmas active in model 2 only
            idx2 = in & (m==2);
            lp(idx2) = lp(idx2) - lambda * theta(idx2,2) - lambda * theta(idx2,4);

            % bounds check on sigma params for model 2 only
            idx2_out_of_bounds = idx2 & ~(theta(:,2) >= limits(2,1) & theta(:,2) <= limits(2,2) ...
                & theta(:,4) >= limits(4,1) & theta(:,4) <= limits(4,2));
            in(idx2_out_of_bounds) = false;

            lp(~in) = -Inf;

        end

        %% distance

        function cmdist = CvMdivergence(sim, obs)
            dists = zeros(2,1);
            for i=1:3
                [~,~,dists(i)] = cmtest2(sim(:,i), obs(:,i));
            end
            cmdist = norm(dists);
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

            activeCols = {[1 3], [1 2 3 4]};
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


        %% prior sampler uniform

        function theta = priorSamplerFn_uniform(N)
            % Draw N parameter vectors from prior
            % abc_smc_adaptive expects theta to be N x D
            % ie each particle is one row; columns are parameters

            D = 4;
            % theta = zeros(N, D);
            % theta will be in this order:
            % 1. muXenter
            % 2. sigmaXenter
            % 3. muXstay
            % 4. sigmaXstay

            limits = uniformPriorLimits();
            limitRange = diff(limits,1,2);
            if size(limits,1) ~= D
                error('number of dimensions in prior limits does not match')
            end
            theta = rand(N,D).*repmat(limitRange',N,1) + repmat(limits(:,1)',N,1);

        end


        %% prior log prob uniform

        function lp = priorLogPdfFn_uniform(theta)
            % Log prior pdf for NxD theta
            %
            % Must return -Inf for impossible values
            [N,D] = size(theta);
            lp = zeros(N,1);
            limits = uniformPriorLimits();

            in = true(N,1);

            for i=1:D
                in = in & theta(:,i) >= limits(i,1) & theta(:,i) <= limits(i,2);
            end

            lp(~in) = -Inf;
            lp(in)  = 0;
        end


    end
end
