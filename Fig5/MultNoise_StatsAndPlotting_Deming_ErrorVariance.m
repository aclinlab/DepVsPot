
%% Statistical Analysis & Plotting of MBON responses
%
% 
% Per-fly analysis script for detecting signal-dependent noise in paired MBON(MBONg2, MBONg5) 
% calcium imaging responses. Run once per fly; results accumulate across
% flies in a shared summary file for population analysis.
%
% Requires DetailedData.mat files for each fly (MultNoise_DataExtraction.m), 
% containing all relevant variables.
%
% Analysis pipeline:
%   1. Deming (TLS) regression on paired MBON responses to extract shared
%      latent drive (proj_sign) and orthogonal residuals (r_orth)
%   2. Two measures of signal-noise scaling:
%       - Spearman rho           : monotonic relationship (rank-based)
%       - OLS slope              : linear scaling magnitude
%       
%   3. Permutation null          : shuffles residuals, preserving stimulus
%                                  means and inter-neuron noise correlation
%   4. Matched shot noise control: TLS pipeline run on size-matched
%                                  split-half photon noise samples
%   5. Population statistics     : Wilcoxon tests, paired permutation tests
%                                  and Monte Carlo null distribution across
%                                  flies
% 
%   Output:
%       - Per-fly results appended to All_Flies_Summary_Deming_WithinVar.mat
%       - Dashboard figure (.fig) with per-fly and population plots saved
%         as SVGs
%       - Timestamped command window log
%
% Dependencies: tls_projection, paired_perm (defined as local functions)
%             
% See also: MultNoise_DataExtraction.m
%
% Written by Katharina Bergmann

%% -------------------------------------------------------------------------
%               1. LOAD FLY DATA
% --------------------------------------------------------------------------

clear; clc; close all;

% Select the shared output folder - must be the SAME folder across all flies
% so that All_Flies_Summary.mat accumulates correctly
fprintf('\nPlease select the "Multiplicative Noise Analysis Results" folder in the dialogue window.\n');
baseFolder = uigetdir(pwd, 'Select output folder');
if baseFolder == 0
    error('User cancelled folder selection');
end
[~,folderName] = fileparts(baseFolder);
if ~strcmp(folderName, 'Multiplicative Noise Analysis Results')
    warning('Selected folder is not named "Multiplicative Noise Analysis Results"');
end

% --- Shared output folder  ---
masterOutputFolder = fullfile(baseFolder, 'Analysis_Output_Deming_WithinVar');
figsFolder         = fullfile(masterOutputFolder, 'Figures_Deming_WithinVar');
diaryFolder        = fullfile(masterOutputFolder, 'Command Log');

if ~exist(masterOutputFolder, 'dir') 
    mkdir(masterOutputFolder); 
end

if ~exist(figsFolder, 'dir')        
    mkdir(figsFolder);         
end

if ~exist(diaryFolder, 'dir')        
    mkdir(diaryFolder);        
end

% --- Load fly data ---
fprintf('\nSelect the "Detailed Data" file for the fly you want to analyse.\n');
[fileName, pathName] = uigetfile('*.mat', 'Select the DetailedData file');
if isequal(fileName, 0); return; end
load(fullfile(pathName, fileName));

% Extract flyID for naming figures and indexing summary
flyID = strrep(fileName, '_DetailedData.mat', '');
fprintf('Loaded fly: %s\n', flyID);

% --- Save timestamped command window log
diaryTimestamp = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
diaryFile      = fullfile(diaryFolder, ['MBON_Analysis_Log_' flyID '_' diaryTimestamp '.txt']);
diary(diaryFile);
diary on;
fprintf('Command window log started: %s\n', diaryFile);

% --- Extract VolumeRate from metadata ---
volumeRates = cellfun(@(SI) SI.hRoiManager.scanVolumeRate, all_SI);

% Sanity check: are they all the same?
if numel(unique(volumeRates)) > 1
    warning('scanVolumeRate differs across trials for this fly — check acquisition settings.');
end

fs = volumeRates(1);   
frequency = [2 4 10 20 40];


%% ------------------------------------------------------------------------- 
%               2. STATISTICAL ANALYSIS - SHOT NOISE AND MBON DATA 
% --------------------------------------------------------------------------
% This part performs Deming regression, a generalisation of total least squares (TLS) regression),
% that accounts for unequal error variances in the two variables, captured by lambda (ratio of mean within-stimulus variances)
% on the two MBONg2/MBONg5 and MBONg2A/B, and MBONg5AB/B responses distributions. 
% Unlike OLS, it minimises orthogonal, NOT vertical distances and is applicable
% when you have noise in both X and Y, treating both variables equally and not assuming
% one is dependent on the other.
% Orthogonal residual = proxy for trial-to-trial variability/noise
% magnitude
% TLS projection = proxy for signal strength /response amplitude

% Relationshop between response amplitude and noise magnitude are tested with
% two measures:
% 1) Spearman correlation coeffcient (rho)- test for monotonicity of relationship
% 2) Linear OLS slope - quantifies linear scaling of noise magnitude with
% response amplitude
% 

% --------------- 12) DEMING REGRESSION -------------------------------------
% Regression axis is used as proxy for latent shared drive
% proj_sign  = signed projection onto TLS axis
% r_orth = absolute orthogonal residuals (euclidean distance)
% r_orth_scaled = corrected for split-half ROIs 

fprintf('\n--- Running Deming regression for this fly  ---\n');

% Flatten 3D shot noise matrices into 1D vectors for TLS
A1_long = reshape(all_shotnoise_A1, [], 1);
B1_long = reshape(all_shotnoise_B1, [], 1);
A2_long = reshape(all_shotnoise_A2, [], 1);
B2_long = reshape(all_shotnoise_B2, [], 1);

% ---------- SHOT NOISE MBON1 (gamma2) ----------

lambda_s1 = 1.0; % Simulated from same data, ratio of within-stimulus variance should be approx. 1
[A1, B1] = deal(A1_long(:), B1_long(:));
[v_tls_s1, proj_sign_s1, r_orth_s1]...
    = tls_projection(A1, B1, 'Shot MBONap', false, lambda_s1);


% ---------- SHOT NOISE MBON2 (gamma5) ----------

lambda_s2 = 1.0; % Simulated from same data, ratio of within-stimulus variance should be approx. 1
[A2, B2] = deal(A2_long(:), B2_long(:));
[v_tls_s2, proj_sign_s2, r_orth_s2] ...
    = tls_projection(A2, B2, 'Shot MBONav', false, lambda_s2);

% ---------  ACTUAL DATA ------------------------

% ---  Calculate noise (error variance) based lambda for Deming regression ---
% Calculate within-stimulus variance (noise only)
% var(..., 0, 2) calculates variance across trials for each stimulus
within_var_ROI1 = mean(var(all_amplitudes_ROI1, 0, 2), 'omitnan');
within_var_ROI2 = mean(var(all_amplitudes_ROI2, 0, 2), 'omitnan');

% Lambda represents the ratio of the noise floor in ROI2 vs ROI1
lambda_act = within_var_ROI2 / within_var_ROI1;

fprintf('Calculated Noise Lambda (ROI2/ROI1): %.3f\n', lambda_act);

amps1 = reshape(all_amplitudes_ROI1, [], 1);
amps2 = reshape(all_amplitudes_ROI2, [], 1);

% TLS projection
[v_tls_actual, proj_sign_actual, r_orth_actual] ...
        = tls_projection(amps1, amps2, 'Actual MBON', false, lambda_act);

% get data range for plotting shot noise on same axes as actual data for
% more intuitive visualisation
buffer_x = range(proj_sign_actual) * 0.05;
buffer_y = range(r_orth_actual) * 0.05;
shared_xlim = [min(proj_sign_actual) - buffer_x, max(proj_sign_actual) + buffer_x];
shared_ylim = [min(r_orth_actual) - buffer_y, max(r_orth_actual) + buffer_y];


% --- Defensive Checks for subsequent statistical analysis ----

% tls_projection must return finite outputs for all inputs
assert(all(isfinite(proj_sign_s1) & isfinite(r_orth_s1)), ...
    'tls_projection returned non-finite values for Shot MBONap');
assert(all(isfinite(proj_sign_s2) & isfinite(r_orth_s2)), ...
    'tls_projection returned non-finite values for Shot MBONav');
assert(all(isfinite(proj_sign_actual) & isfinite(r_orth_actual)), ...
    'tls_projection returned non-finite values for Actual MBON');



% --------------- 13) STATISTICAL ANALYSIS OF SIGNAL-NOISE RELATIONSHIP ---
%
%         1. Spearman correlation - rho
%            Monotonic relationship between response amplitude and noise magnitude
%
%         2. Linear OLS regression - slope
%            OLS slope of response amplitude vs noise magnitude quantifies how strongly noise scales
%            with signal. If actual slope > shot noise slope, this suggests biological
%            signal-dependent noise on top of photon shot noise


% --- Scale shot noise residuals to full-ROI equivalent ---
% 
% Relative shot noise (as a proportion of signal) scales as 1/√N_photons, 
% and N_photons ∝ N_pixels (assuming roughly uniform brightness across
% ROI). Therefore:

% a) A half-ROI is noisier than full ROI by √2:
%   sigma_half / sigma_full = (1/√(N/2)) / (1/√N) = √N / √(N/2) = √2 
%
% 
% b) The halves share same underlying signal, but have indepedent photon
% noise, so signal cancels:
%   A = S + e_A, B = S + e_B
%   A - B = e_A - e_B
% Variances add for independent noise:
%   Var(A - B) = sigma_A^2 + sigma_B^2 = 2 * sigma_half^2
%   SD(A - B) = √2 x sigma_half
% 
% c) Both halves measure same signal with equal noise, lambda = 1 and TLS
% axis lies on 45-degree line
%   v_tls = [1/√2 ; 1/√2]
%   v_perp = [-1/√2 ; 1/√2]
%   Projecting onto v_perp gives
%   r_orth = |A - B| / √2  
% 
% d) Combining b) and c):
%   SD(r_orth) = SD(A − B) / √2 = sigma_half
%
% e) Substituting a)
%   SD(r_orth) = √2 * sigma_full
%  
% To make shot noise comparable to actual data (full ROI), divide every
% r_orth by √2, which puts SD(r_orth) on full ROI scale.
% 
% n_pix is saved in loaded DetailedData file 
%
% NB: Scaling applied only for OLS, Spearman rho uses
% raw r_orth (rank-based measure, scaling has no effect)

n_pix_full_s1 = n_pix_1;
n_pix_half_s1 = floor(n_pix_full_s1/2);
r_orth_s1_scaled = r_orth_s1 / sqrt(n_pix_full_s1 / n_pix_half_s1); % scale = sqrt2

n_pix_full_s2 = n_pix_2;
n_pix_half_s2 = floor(n_pix_full_s2/2);
r_orth_s2_scaled = r_orth_s2 / sqrt(n_pix_full_s2 / n_pix_half_s2); % scale = sqrt2


% --- a) SPEARMAN CORRELATION ---
% Is there a monotonic relationship between response amplitude and noise
% magnitude? 
% NB: Rank-based, so no scaling of shot noise data needed

fprintf('\n--- Statistical Analysis Results ---\n');

% SHOT NOISE MBON1
[rho_shot1, p_shot1] = corr(proj_sign_s1, r_orth_s1, ...
    'Type','Spearman');
fprintf('Shot MBONap: Spearman rho = %.3f, p = %.4g, n = %d\n', ...
    rho_shot1, p_shot1, numel(proj_sign_s1));

% SHOT NOISE MBON2 
[rho_shot2, p_shot2] = corr(proj_sign_s2, r_orth_s2, ...
    'Type','Spearman');
fprintf('Shot MBONav: Spearman rho = %.3f, p = %.4g, n = %d\n', ...
    rho_shot2, p_shot2, numel(proj_sign_s2));

%  ACTUAL MBON DATA 
[rho_actual, p_actual] = corr(proj_sign_actual, r_orth_actual, ...
        'Type','Spearman');
fprintf('Actual MBON: Spearman rho = %.3f, p = %.4g, n = %d\n', ...
        rho_actual, p_actual, numel(proj_sign_actual));



% --- b) LINEAR OLS ---
% Quantify scaling of noise magnitude with response amplitude
% X (predictor) → proj_sign_ → response amplitude
% Y (response) → r_orth_ → noise magnitude

% SHOT NOISE 1
mdl_shot1 = fitlm(proj_sign_s1, r_orth_s1_scaled);
slope_shot1 = mdl_shot1.Coefficients.Estimate(2);

% SHOT NOISE 2
mdl_shot2 = fitlm(proj_sign_s2, r_orth_s2_scaled);
slope_shot2 = mdl_shot2.Coefficients.Estimate(2);

% ACTUAL DATA
mdl_act = fitlm(proj_sign_actual, r_orth_actual);
slope_act = mdl_act.Coefficients.Estimate(2);

% Print ouput
fprintf('OLS slopes (linear): Shot1=%.3f | Shot2=%.3f | Actual=%.3f\n', ...
    slope_shot1, slope_shot2, slope_act);


% --------------- 14a) PERMUTATION CONTROL - PER-FLY LEVEL ----------------------------
% 
% Could the observed coupling between response amplitude and noise magnitude arise 
% from analysis pipeline alone?
%
% This control breaks the stimulus-noise coupling by shuffling 2D residual vectors across
% stimulus levels while preserving:
%   - Empirical noise distribution (real residuals, not parametric)
%   - Inter-neuron noise correlation (residual pairs shuffled jointly)
%   - Axis-estimation artifact (full TLS pipeline run on each permutation)
%
% Null hypothesis: residual magnitude is exchangeable across stimulus
% levels, i.e., noise does not scale with response amplitude
% Logic:
%   1. Compute per-stimulus mean response for each neuron
%   2. Compute residuals (data - stimulus mean) per trial per stimulus
%   3. Pool all 50 2D residual vectors
%   4. Randomly reassign residuals to stimulus means → null dataset
%   5. Re-estimate lambda from the null dataset, fit TLS, compute statistics
%   6. Store null distributions of Spearman rho and OLS slope
%   7. Empirical one-sided p-values against each fly's null distribution 


nPerm_null  = 10000;

fprintf('\n--- Running permutation null control (n=%d) ---\n', nPerm_null);


% --- Dynamic sizing to handle varying trial numbers ---
[nStim, nTrials] = size(all_amplitudes_ROI1);   % Automatically detect [5, 9] or [5, 10]
nTotal           = numel(all_amplitudes_ROI1);  % Total points (e.g., 45 or 50)

% --- Reshape actual amplitudes into [nStim x nTrials] ---
amp1_mat = all_amplitudes_ROI1;   % [5 x 10] stim x trials
amp2_mat = all_amplitudes_ROI2;   % [5 x 10] stim x trials

% --- Step 1: Per-stimulus means [nStim x 1] ---
mean1_per_stim = mean(amp1_mat, 2);   % [5 x 1]
mean2_per_stim = mean(amp2_mat, 2);   % [5 x 1]

% --- Step 2: Residuals per trial per stimulus ---
resid1_mat = amp1_mat - mean1_per_stim;   % [5 x 10]
resid2_mat = amp2_mat - mean2_per_stim;   % [5 x 10]

% --- Step 3: Pool all 50 2D residual vectors ---
resid1_pool   = resid1_mat(:);              % [50 x 1] Trial 1 s1,s2...,s5, Trial 2 s1,s2...
resid2_pool   = resid2_mat(:);              % [50 x 1]
resid_pool_2D = [resid1_pool, resid2_pool]; % [50 x 2] paired 2D residuals

% --- Replicate stimulus means for all trials ---
% repmat replicates stimulus means nTrials times to match the column-major
% flattening of resid_mat(:), where all stimuli for trial 1 come first,
% then all stimuli for trials 2, etc
means1_rep = repmat(mean1_per_stim, nTrials, 1);   % [50 x 1] Same structure as resid1_pool: Trial 1 s1,s2...,s5, Trial 2 s1,s2...
means2_rep = repmat(mean2_per_stim, nTrials, 1);   % [50 x 1]

% --- Preallocate null distributions ---
null_rho_actual   = nan(nPerm_null, 1);
null_slope_actual = nan(nPerm_null, 1);


% --- Step 4-5: Permutation loop ---
for p = 1:nPerm_null

    % Randomly shuffle 2D residual vectors jointly (preserves inter-neuron pairs)
    shuf_idx      = randperm(nTotal);
    resid_shuf_2D = resid_pool_2D(shuf_idx, :);   % [50 x 2]

    % Reconstruct null dataset: stimulus mean + shuffled residual
    null_amp1 = means1_rep + resid_shuf_2D(:, 1);   % [50 x 1]
    null_amp2 = means2_rep + resid_shuf_2D(:, 2);   % [50 x 1]

    % --- Calculate null_lambda for each permutation ---
    % For strict exchangeability should not re-use lambda_act 
    % Shuffled residuals change within-stimulus variance structure, so noise ratio should be
    % re-estimated for each null dataset.
    % NB: Test first whether this results in degenerate lambda values that destabilise TLS axis
    % - if yes, just use lambda_act (see below)

    % Reshape null amplitudes back to [nStim x nTrials] to compute within-stimulus variance
    null_amp1_mat = reshape(null_amp1, nStim, nTrials);
    null_amp2_mat = reshape(null_amp2, nStim, nTrials);

    % Within-stimulus variance for this permutation
    null_within_var1 = mean(var(null_amp1_mat, 0, 2), 'omitnan');
    null_within_var2 = mean(var(null_amp2_mat, 0, 2), 'omitnan');

    % Calculate null_lambda
    null_lambda = null_within_var2 / null_within_var1;

    % Run full TLS pipeline on null dataset 
    [~, null_proj, null_rorth] = tls_projection(null_amp1, null_amp2, '', true, null_lambda);
    
    % Permutation test should have similar points as real data, but still possible
    % to reject degenerate permutations. If get warning a lot -> null_lambda
    % is regularly unstable, better use lambda_actual
    valid_null = isfinite(null_proj) & isfinite(null_rorth);
    min_valid = round(0.9 * nTotal);   % require 90% of points to be finite
    
    if sum(valid_null) < min_valid
        warning('Permutation %d: only %d/%d finite values (null_lambda=%.4f). Skipping.', ...
        p, sum(valid_null), nTotal, null_lambda);
      continue
    end

    % Passed 90% threshold, now extract valid points
    np_clean = null_proj(valid_null); 
    nr_clean = null_rorth(valid_null);  

    % Spearman rho
    null_rho_actual(p) = corr(np_clean, nr_clean, ...
                              'Type', 'Spearman');

    % OLS slope 
    coeffs_slope = polyfit(np_clean, nr_clean, 1);
    null_slope_actual(p) = coeffs_slope(1);

end

% --- Step 6: Extract finite permutation values for both statistics ---
null_rho   = null_rho_actual(isfinite(null_rho_actual));
null_slope = null_slope_actual(isfinite(null_slope_actual));

fprintf('Valid permutations: rho=%d | slope=%d (of %d)\n', ...
    numel(null_rho), numel(null_slope), nPerm_null);


fprintf('Observed actual slope:   %.4f\n', slope_act);
fprintf('Permutation null slope median: %.4f\n', median(null_slope));
 

% ---  Step 7: Empirical p-values (one-sided, add-one convention) ---
n_exc_rho   = sum(null_rho   >= rho_actual);
n_exc_slope = sum(null_slope >= slope_act);

p_perm_rho_actual   = (1 + n_exc_rho)   / (1 + numel(null_rho));
p_perm_slope_actual = (1 + n_exc_slope) / (1 + numel(null_slope));

fprintf('\nPermutation null results:\n');
fprintf('  Rho:   observed=%.3f | null median=%.3f | %d/%d null >= obs | p=%.2g\n', ...
    rho_actual, median(null_rho),   n_exc_rho,   numel(null_rho),   p_perm_rho_actual);
fprintf('  Slope: observed=%.3f | null median=%.3f | %d/%d null >= obs | p=%.2g\n', ...
    slope_act,  median(null_slope), n_exc_slope, numel(null_slope), p_perm_slope_actual);

% NB: report medians here (not mean as in population analysis) as per fly
% nulls are skewed, compared to relatively symmetric population null

% BUT print means in summary violin plot, so its consistent within figure


% --------------- 14b) SAMPLE-SIZE MATCHED SHOT NOISE CONTROL - PER-FLY LEVEL -------------
% Logic:
%   Actual data has 50 points (5 stimuli x 10 trials). The full shot noise
%   control pools all 5 x 1000 x 10 = 50,000 simulated points into a single fit, 
%   so its TLS axis is estimated far more stably than the axis fitted to 50 real points.
%   A steeper OLS slope in actual data could therefore
%   reflect small-sample estimation volatility rather than biology.
%
%   This control removes that asymmetry by refitting the shot noise null at
%   the actual sample size: each pre-generated shot noise realization is
%   treated as its own 50-point dataset (5 stimuli x 10 trials = one iteration slice),
%   undergoes full TLS + OLS pipeline independently, and builds a slope distribution under the shot noise
%   null at matched sample size.
%
% Method - Monte Carlo null:

% Shot noise array structure: [5 x 1000 x 10]
%   dim 1 = stimuli (5)
%   dim 2 = iterations (1000)
%   dim 3 = trials (10)
%
%   The loop walks deterministically over k = 1:1000, taking slice
%   (:, k, :) -> [5 x 10] -> reshape to [50 x 1], mirroring actual data
%   structure. Every realization is used exactly once, in order -> creates a Monte
%   Carlo null over 1000 independent realizations of the shot noise
%   generative model.
%   
%   Lambda is recomputed within each slice from its own 50 points, exactly
%   as it is for the actual data (unlike for pooled shot noise data, where lambda set to 1)
%
% ROI combination:
%   ROI1 and ROI2 are run separately (preserving per-ROI lambda / pixel
%   scaling), then their slope distributions are concatenated into one combined
%   null. NB: this makes the pooled null a mixture of two
%   per-ROI nulls rather than 2000 draws of a single quantity; p-values
%   below are therefore against "shot noise in a randomly chosen ROI".
%
% Output per fly (saved in section 15 alongside all other per-fly fields):
%   - slope_dist_sn_matched_pooled_all : pooled OLS slope dist (cell, 2000 x 1 per fly)
%   - slope_sn_matched_median_all      : median of pooled slope distribution
%   - p_sn_matched_slope_all           : empirical, one-sided p-value
%
%   - rho_dist_sn_matched_pooled_all   : " ------------------------"
%   - rho_sn_matched_median_all
%   - p_sn_matched_rho_all
%
% Population test (also saved in section 15):
%   Equivalent to population null for permutation control
%   sample one value per fly from pooled distribution, take mean, repeat
%   100,000 times → null distribution of population mean slope.

nIter_matched = 1000;   % number of iteration slices (dim 2) in shot noise arrays

% --- Confirm array dimensionality ---
[nStim_sn, nIter_sn_dim, nTrials_sn] = size(all_shotnoise_A1);

assert(nStim_sn   == 5,                     'Expected 5 stimuli in dim 1 of all_shotnoise_A1');
assert(nTrials_sn == 10 || nTrials_sn ==9,  'Expected 9/10 trials in dim 3 of all_shotnoise_A1');
assert(nIter_sn_dim >= nIter_matched,       'Requested %d iterations but only %d available', ...
                                       nIter_matched, nIter_sn_dim);

fprintf('\n--- Sample-size matched shot noise control (n=%d x %d-pt samples per ROI) ---\n',...
    nIter_matched, nStim_sn * nTrials_sn);

% --- Pre-allocate per-ROI distributions ---
slope_dist_sn_matched_1 = nan(nIter_matched, 1);
slope_dist_sn_matched_2 = nan(nIter_matched, 1);
rho_dist_sn_matched_1 = nan(nIter_matched, 1);
rho_dist_sn_matched_2 = nan(nIter_matched, 1);

min_pts_sn = 45;     % require >= 45/50 finite points (90%) per iteration

% Scale factors to convert half-ROI residuals to full-ROI equivalent.
% Same correction as applied to the pooled shot noise control in section 13
% (see derivation there). proj is unaffected - both halves carry the same
% signal - so only r_orth needs scaling.
scale_s1 = sqrt(n_pix_full_s1 / n_pix_half_s1);   % = sqrt(2)
scale_s2 = sqrt(n_pix_full_s2 / n_pix_half_s2);   % = sqrt(2)

for k = 1:nIter_matched

    % ---- ROI 1 ----
    sn1_A = reshape(all_shotnoise_A1(:, k, :), [], 1);   % [50 x 1]
    sn1_B = reshape(all_shotnoise_B1(:, k, :), [], 1);

    sn1_A_mat    = reshape(sn1_A, nStim_sn, nTrials_sn);
    sn1_B_mat    = reshape(sn1_B, nStim_sn, nTrials_sn);
    wv_sn1_A     = mean(var(sn1_A_mat, 0, 2), 'omitnan');
    wv_sn1_B     = mean(var(sn1_B_mat, 0, 2), 'omitnan');
    lambda_sn1_k = wv_sn1_B / wv_sn1_A;

    [~, proj_sn1, r_sn1] = tls_projection(sn1_A, sn1_B, '', true, lambda_sn1_k);

    valid_sn1 = isfinite(proj_sn1) & isfinite(r_sn1);
    if sum(valid_sn1) >= min_pts_sn
        p1 = proj_sn1(valid_sn1);
        r1 = r_sn1(valid_sn1);
        r1_scaled = r1 / scale_s1; % full ROI equivalent

        c1 = polyfit(p1, r1_scaled, 1);
        slope_dist_sn_matched_1(k) = c1(1);

        % NB: rho is rank-based, so no scaling needed
        rho_dist_sn_matched_1(k) = corr(p1, r1, 'Type', 'Spearman');
    end

    % ---- ROI 2 ----
    sn2_A = reshape(all_shotnoise_A2(:, k, :), [], 1);
    sn2_B = reshape(all_shotnoise_B2(:, k, :), [], 1);

    sn2_A_mat    = reshape(sn2_A, nStim_sn, nTrials_sn);
    sn2_B_mat    = reshape(sn2_B, nStim_sn, nTrials_sn);
    wv_sn2_A     = mean(var(sn2_A_mat, 0, 2), 'omitnan');
    wv_sn2_B     = mean(var(sn2_B_mat, 0, 2), 'omitnan');
    lambda_sn2_k = wv_sn2_B / wv_sn2_A;

    [~, proj_sn2, r_sn2] = tls_projection(sn2_A, sn2_B, '', true, lambda_sn2_k);

    valid_sn2 = isfinite(proj_sn2) & isfinite(r_sn2);
    if sum(valid_sn2) >= min_pts_sn
        p2 = proj_sn2(valid_sn2);
        r2 = r_sn2(valid_sn2);
        r2_scaled = r2 / scale_s2; % full ROI equivalent

        c2 = polyfit(p2, r2_scaled, 1);
        slope_dist_sn_matched_2(k) = c2(1);

        rho_dist_sn_matched_2(k) = corr(p2, r2, 'Type', 'Spearman');
    end

end

% --- Pool ROI1 + ROI2 into single combined null distribution ---
slope_sn_matched_1_finite = slope_dist_sn_matched_1(isfinite(slope_dist_sn_matched_1));
slope_sn_matched_2_finite = slope_dist_sn_matched_2(isfinite(slope_dist_sn_matched_2));
rho_sn_matched_1_finite = rho_dist_sn_matched_1(isfinite(rho_dist_sn_matched_1));
rho_sn_matched_2_finite = rho_dist_sn_matched_2(isfinite(rho_dist_sn_matched_2));

slope_sn_matched_pooled = [slope_sn_matched_1_finite; slope_sn_matched_2_finite];
rho_sn_matched_pooled = [rho_sn_matched_1_finite; rho_sn_matched_2_finite];

fprintf('Valid iterations: ROI1=%d | ROI2=%d | pooled=%d (of %d requested)\n', ...
    numel(slope_sn_matched_1_finite), numel(slope_sn_matched_2_finite), ...
    numel(slope_sn_matched_pooled),   nIter_matched);

% --- Per fly summary statistics ---
slope_sn_matched_median = median(slope_sn_matched_pooled);
rho_sn_matched_median = median(rho_sn_matched_pooled);

% --- Empirical p-values (one-sided, add-one convention) ---
n_exc_sn_slope = sum(slope_sn_matched_pooled >= slope_act);
n_exc_sn_rho = sum(rho_sn_matched_pooled >= rho_actual);

p_sn_matched_slope = (1 + n_exc_sn_slope) / (1 + numel(slope_sn_matched_pooled));
p_sn_matched_rho = (1 + n_exc_sn_rho) / (1 + numel(rho_sn_matched_pooled));

fprintf('\nMatched-sample shot noise control results:\n');
fprintf('  Slope: pooled median=%.4f | actual=%.4f | %d/%d null >= obs | p=%.2g\n', ...
    slope_sn_matched_median, slope_act, n_exc_sn_slope, numel(slope_sn_matched_pooled), p_sn_matched_slope);
fprintf('  Rho: pooled median=%.4f | actual=%.4f | %d/%d null >= obs | p=%.2g\n', ...
    rho_sn_matched_median, rho_actual, n_exc_sn_rho, numel(rho_sn_matched_pooled), p_sn_matched_rho);

% NB: report medians here (not mean as in population analysis) as per fly
% nulls are skewed, compared to relatively symmetric population null



% -------------- 15) SAVE PER-FLY ACCUMULATORS IN MASTER FILE --------

summaryFile = fullfile(masterOutputFolder, 'All_Flies_Summary_Deming_WithinVar.mat');

% Load existing summary or create fresh struct
if exist(summaryFile, 'file')
    S = load(summaryFile);
    fprintf('Loaded summary: %d flies: %s\n', numel(S.fly_IDs_all), strjoin(S.fly_IDs_all, ', '));
    if isfield(S, 'perm_null_slope_dist')
        fprintf('perm_null_slope_dist length: %d\n', numel(S.perm_null_slope_dist));
        fprintf('Non-empty entries: %d\n', sum(cellfun(@(x) ~isempty(x), S.perm_null_slope_dist)));
    end
else
    S = struct();
    S.fly_IDs_all = {};
    fprintf('No existing summary found. Starting fresh.\n');
end

% --- Find or create index for this fly ---
idx = find(strcmp(S.fly_IDs_all, flyID), 1);
if isempty(idx)
    idx = numel(S.fly_IDs_all) + 1;
    S.fly_IDs_all{idx} = flyID;
else
    warning('Fly %s already exists in summary. Overwriting previous entry.', flyID);
end

% --- Initialise all scalar accumulator fields if missing or too short ---
fields = {'rho_shot1_all','rho_shot2_all','rho_actual_all', ...
          'p_shot1_all','p_shot2_all','p_actual_all', ...
          'slope_shot1_all','slope_shot2_all','slope_actual_all', ...
          'perm_null_rho_actual', 'perm_null_slope_actual', ...
          'perm_null_rho_median', 'perm_null_slope_median', ...
          'slope_sn_matched_median_all', 'rho_sn_matched_median_all',...
          'p_sn_matched_slope_all', 'p_sn_matched_rho_all'};

for f = fields
    if ~isfield(S, f{1})
        S.(f{1}) = nan(1, idx);
    elseif numel(S.(f{1})) < idx
        S.(f{1})(end+1:idx) = NaN;
    end
end


% --- Initialise cell array fields ---
% Uses cell(1, idx) to guarantee all positions 1:idx are properly 
% initialised as [] rather than leaving gaps from direct assignment (NB: numeric array explicitly fills gaps with NaNs)
cell_fields = {'perm_null_rho_dist', 'perm_null_slope_dist', ...
    'proj_actual_dist', 'r_orth_actual_dist', 'slope_dist_sn_matched_pooled_all', ...
    'rho_dist_sn_matched_pooled_all', 'amplitudes_ROI1_all', 'amplitudes_ROI2_all'};

for f = cell_fields
    if ~isfield(S, f{1})
        S.(f{1}) = cell(1, idx);      % all positions initialised as []
    elseif numel(S.(f{1})) < idx
        S.(f{1}){idx} = [];           % MATLAB fills gaps as [], not uninitialised
    end
end


% --- Assign this fly's values --- Guards removed as always start afresh
%                                  for new analysis

S.rho_shot1_all(idx)                    = rho_shot1;
S.rho_shot2_all(idx)                    = rho_shot2;
S.rho_actual_all(idx)                   = rho_actual;

S.p_shot1_all(idx)                      = p_shot1;
S.p_shot2_all(idx)                      = p_shot2;
S.p_actual_all(idx)                     = p_actual;

S.slope_shot1_all(idx)                  = slope_shot1;
S.slope_shot2_all(idx)                  = slope_shot2;
S.slope_actual_all(idx)                 = slope_act;

S.perm_null_rho_actual(idx)             = p_perm_rho_actual;
S.perm_null_slope_actual(idx)           = p_perm_slope_actual;
S.perm_null_rho_median(idx)             = median(null_rho);
S.perm_null_slope_median(idx)           = median(null_slope);

S.perm_null_rho_dist{idx}               = null_rho;
S.perm_null_slope_dist{idx}             = null_slope;

S.proj_actual_dist{idx}                 = proj_sign_actual;
S.r_orth_actual_dist{idx}               = r_orth_actual;

S.slope_sn_matched_median_all(idx)      = slope_sn_matched_median;
S.rho_sn_matched_median_all(idx)        = rho_sn_matched_median;
S.p_sn_matched_slope_all(idx)           = p_sn_matched_slope;
S.p_sn_matched_rho_all(idx)             = p_sn_matched_rho;
S.slope_dist_sn_matched_pooled_all{idx} = slope_sn_matched_pooled;
S.rho_dist_sn_matched_pooled_all{idx}   = rho_sn_matched_pooled;

S.amplitudes_ROI1_all{idx}              = all_amplitudes_ROI1;
S.amplitudes_ROI2_all{idx}              = all_amplitudes_ROI2;


% --- Save per-fly data ---
save(summaryFile, '-struct', 'S');
fprintf('Per-fly data saved to: %s\n', summaryFile);



% ---------------- 16) POPULATION STATISTICS -------------------

datasets = {'Shot MBONap','Shot MBONav','Actual'};
rho_sets = {S.rho_shot1_all, S.rho_shot2_all, S.rho_actual_all};

% --- One-sided Wilcoxon test (against 0) for rho ---

S.pop.rho_wilcoxon = struct([]);

for d = 1:numel(datasets)
    r = rho_sets{d};
    r = r(isfinite(r));

    S.pop.rho_wilcoxon(d).dataset = datasets{d};
    S.pop.rho_wilcoxon(d).n = numel(r);

    if numel(r) >= 6 % smallest possible p = 1/2^n
        [p,~,~] = signrank(r,0,'tail','right');
        S.pop.rho_wilcoxon(d).median_rho = median(r);
        S.pop.rho_wilcoxon(d).p_value = p;
    else
        S.pop.rho_wilcoxon(d).median_rho = NaN;
        S.pop.rho_wilcoxon(d).p_value = NaN;
        warning('Not enough data for Wilcoxon test on %s (n=%d).', datasets{d}, numel(r));
    end
end


% --- Paired Wilcoxon test for slopes ---
% Actual vs Shot1
if numel(S.slope_actual_all) == numel(S.slope_shot1_all) && ...
   numel(S.slope_actual_all) >= 6

    [p1,~,~] = signrank(S.slope_actual_all, S.slope_shot1_all, 'tail','right');
    S.pop.slope_wilcoxon(1).comparison  = 'Actual vs Shot1';
    S.pop.slope_wilcoxon(1).median_diff = median(S.slope_actual_all - S.slope_shot1_all, 'omitnan');
    S.pop.slope_wilcoxon(1).p_value     = p1;
else
    S.pop.slope_wilcoxon(1).comparison  = 'Actual vs Shot1';
    S.pop.slope_wilcoxon(1).median_diff = median(S.slope_actual_all - S.slope_shot1_all, 'omitnan');
    S.pop.slope_wilcoxon(1).p_value     = NaN;
    warning('Too few paired samples for Wilcoxon test: Actual vs Shot1.');
end

% Actual vs Shot2
if numel(S.slope_actual_all) == numel(S.slope_shot2_all) && ...
   numel(S.slope_actual_all) >= 6


    [p2,~,~] = signrank(S.slope_actual_all, S.slope_shot2_all, 'tail','right');
    S.pop.slope_wilcoxon(2).comparison  = 'Actual vs Shot2';
    S.pop.slope_wilcoxon(2).median_diff = median(S.slope_actual_all - S.slope_shot2_all, 'omitnan');
    S.pop.slope_wilcoxon(2).p_value     = p2;
else
    S.pop.slope_wilcoxon(2).comparison  = 'Actual vs Shot2';
    S.pop.slope_wilcoxon(2).median_diff = median(S.slope_actual_all - S.slope_shot2_all, 'omitnan');
    S.pop.slope_wilcoxon(2).p_value     = NaN;
    warning('Too few paired samples for Wilcoxon test: Actual vs Shot2.');
end

% --- Paired permutation tests for slopes ---
% NB empirical null hypothesis - builds null directly from data (random label flips)

nPerm = 100000;
statFun = @mean;

minPairs = 6; % minimum number of paired samples for permutation

% Slope: Actual vs Shot1
if numel(S.slope_actual_all) == numel(S.slope_shot1_all) && ...   % protect against mismatched arrays being passed to paired_perm
   numel(S.slope_actual_all) >= minPairs
    S.pop.slope_paired_perm.p_shot1 = paired_perm( ...
    S.slope_actual_all, S.slope_shot1_all, statFun, nPerm);
else
    S.pop.slope_paired_perm.p_shot1 = NaN;
    warning('Too few paired samples for permutation test: Actual vs Shot1.');
end

% Slope: Actual vs Shot2
if numel(S.slope_actual_all) == numel(S.slope_shot2_all) && ...
   numel(S.slope_actual_all)>= minPairs
    S.pop.slope_paired_perm.p_shot2 = paired_perm( ...
    S.slope_actual_all, S.slope_shot2_all, statFun, nPerm);
else
    S.pop.slope_paired_perm.p_shot2 = NaN;
    warning('Too few paired samples for permutation test: Actual vs Shot2.');
end

S.pop.slope_paired_perm.mean_diff_shot1 = mean(S.slope_actual_all - S.slope_shot1_all,'omitnan');
S.pop.slope_paired_perm.mean_diff_shot2 = mean(S.slope_actual_all - S.slope_shot2_all,'omitnan');

% --- Paired Wilcoxon test for rho: actual vs shot noise ---
if numel(S.rho_actual_all) == numel(S.rho_shot1_all) && ...
   numel(S.rho_actual_all) >= 6

    [p1,~,~] = signrank(S.rho_actual_all, S.rho_shot1_all, 'tail','right');
    S.pop.rho_wilcoxon_paired(1).comparison  = 'Actual vs Shot1';
    S.pop.rho_wilcoxon_paired(1).median_diff = median(S.rho_actual_all - S.rho_shot1_all, 'omitnan');
    S.pop.rho_wilcoxon_paired(1).p_value     = p1;
else
    S.pop.rho_wilcoxon_paired(1).comparison  = 'Actual vs Shot1';
    S.pop.rho_wilcoxon_paired(1).median_diff = median(S.rho_actual_all - S.rho_shot1_all, 'omitnan');
    S.pop.rho_wilcoxon_paired(1).p_value     = NaN;
    warning('Too few paired samples for Wilcoxon test: rho Actual vs Shot1.');
end

if numel(S.rho_actual_all) == numel(S.rho_shot2_all) && ...
   numel(S.rho_actual_all) >= 6

    [p2,~,~] = signrank(S.rho_actual_all, S.rho_shot2_all, 'tail','right');
    S.pop.rho_wilcoxon_paired(2).comparison  = 'Actual vs Shot2';
    S.pop.rho_wilcoxon_paired(2).median_diff = median(S.rho_actual_all - S.rho_shot2_all, 'omitnan');
    S.pop.rho_wilcoxon_paired(2).p_value     = p2;
else
    S.pop.rho_wilcoxon_paired(2).comparison  = 'Actual vs Shot2';
    S.pop.rho_wilcoxon_paired(2).median_diff = median(S.rho_actual_all - S.rho_shot2_all, 'omitnan');
    S.pop.rho_wilcoxon_paired(2).p_value     = NaN;
    warning('Too few paired samples for Wilcoxon test: rho Actual vs Shot2.');
end

% --- Paired permutation test for rho ---
if numel(S.rho_actual_all) == numel(S.rho_shot1_all) && ...
   numel(S.rho_actual_all) >= minPairs
    S.pop.rho_paired_perm.p_shot1 = paired_perm( ...
        S.rho_actual_all, S.rho_shot1_all, statFun, nPerm);
else
    S.pop.rho_paired_perm.p_shot1 = NaN;
    warning('Too few paired samples for permutation test: rho Actual vs Shot1.');
end

if numel(S.rho_actual_all) == numel(S.rho_shot2_all) && ...
   numel(S.rho_actual_all) >= minPairs
    S.pop.rho_paired_perm.p_shot2 = paired_perm( ...
        S.rho_actual_all, S.rho_shot2_all, statFun, nPerm);
else
    S.pop.rho_paired_perm.p_shot2 = NaN;
    warning('Too few paired samples for permutation test: rho Actual vs Shot2.');
end

S.pop.rho_paired_perm.mean_diff_shot1 = mean(S.rho_actual_all - S.rho_shot1_all, 'omitnan');
S.pop.rho_paired_perm.mean_diff_shot2 = mean(S.rho_actual_all - S.rho_shot2_all, 'omitnan');


% --- Population null for permutation control (14a) ---
% Monte Carlo approximation of null distribution of population mean:
% On each draw: sample one value (with replacement) from each fly's permutation null,
% then average across flies. Repeat 100,000 (nDraw_pop) times. 
% Tests whether the observed mean across flies exceeds this null (one-sided).
% For rho and slope

nDraw_pop = 100000;

nFlies = numel(S.fly_IDs_all);

% Check all flies have null distributions stored
has_null_slope = cellfun(@(x) ~isempty(x), S.perm_null_slope_dist);
has_null_rho = cellfun(@(x) ~isempty(x), S.perm_null_rho_dist);
has_slope = isfinite(S.slope_actual_all);
has_rho = isfinite(S.rho_actual_all);

valid_pop_slope = has_null_slope & has_slope; % if not checking both, could silently misalign if fly has one but NaN for the other
valid_pop_rho   = has_null_rho   & has_rho;


if sum(valid_pop_slope) < nFlies
    warning('%d of %d flies missing null slope distribution. Population test uses available flies only.', ...
        nFlies - sum(valid_pop_slope), nFlies);
end

fly_nulls_slope = S.perm_null_slope_dist(valid_pop_slope);   % cell array, one vector per fly
fly_nulls_rho   = S.perm_null_rho_dist(valid_pop_rho);

nFlies_valid = numel(fly_nulls_slope);

% Observed population mean 
obs_mean_slope = mean(S.slope_actual_all(valid_pop_slope), 'omitnan');
obs_mean_rho   = mean(S.rho_actual_all(valid_pop_rho),     'omitnan');

% Monte Carlo null distributions of the population mean slope and rho
pop_null_mean_slopes = nan(nDraw_pop, 1);
pop_null_mean_rhos   = nan(nDraw_pop, 1);

for b = 1:nDraw_pop
    % Sample one null value per fly, then average
    pop_null_mean_slopes(b) = mean(cellfun(@(v) v(randi(numel(v))), fly_nulls_slope));
    pop_null_mean_rhos(b)   = mean(cellfun(@(v) v(randi(numel(v))), fly_nulls_rho));
end

n_exc_pop_slope = sum(pop_null_mean_slopes >= obs_mean_slope);
n_exc_pop_rho   = sum(pop_null_mean_rhos   >= obs_mean_rho);

p_pop_slope = (1 + n_exc_pop_slope) / (1 + nDraw_pop);
p_pop_rho   = (1 + n_exc_pop_rho)   / (1 + nDraw_pop);

% 95th percentile threshold
threshold_95 = prctile(pop_null_mean_slopes, 95);

fprintf('\n--- Population permutation test ---\n');
fprintf('  Slope: N=%d | observed=%.4f | null mean=%.4f | %d/%d null >= obs | p=%.2g\n', ...
    numel(fly_nulls_slope), obs_mean_slope, mean(pop_null_mean_slopes), ...
    n_exc_pop_slope, nDraw_pop, p_pop_slope);
fprintf('  Rho:   N=%d | observed=%.4f | null mean=%.4f | %d/%d null >= obs | p=%.2g\n', ...
    numel(fly_nulls_rho), obs_mean_rho, mean(pop_null_mean_rhos), ...
    n_exc_pop_rho, nDraw_pop, p_pop_rho);

% Save population results
S.pop.slope_perm.obs_mean             = obs_mean_slope;
S.pop.slope_perm.p_pop_slope          = p_pop_slope;
S.pop.slope_perm.threshold_95         = threshold_95;
S.pop.slope_perm.null_median          = median(pop_null_mean_slopes);
S.pop.slope_perm.null_mean            = mean(pop_null_mean_slopes);
S.pop.slope_perm.nFlies               = nFlies_valid;
S.pop.slope_perm.null_dist            = pop_null_mean_slopes;
S.pop.slope_perm.n_exceed             = n_exc_pop_slope;
S.pop.slope_perm.nDraw                = nDraw_pop;

S.pop.rho_perm.obs_mean               = obs_mean_rho;
S.pop.rho_perm.p_pop_rho              = p_pop_rho;
S.pop.rho_perm.threshold_95           = prctile(pop_null_mean_rhos, 95);
S.pop.rho_perm.n_exceed               = n_exc_pop_rho;
S.pop.rho_perm.nDraw                  = nDraw_pop;
S.pop.rho_perm.null_median            = median(pop_null_mean_rhos);
S.pop.rho_perm.null_mean              = mean(pop_null_mean_rhos);
S.pop.rho_perm.nFlies                 = numel(fly_nulls_rho);
S.pop.rho_perm.null_dist              = pop_null_mean_rhos;



% --- Population null for sample-size matched shot noise control (14b) ---
% Same Monte Carlo procedure as 14a, but each fly's null is its simulated 
% shot-noise slope distribution (fitted at actual sample size) rather than
% its permutation null.
% Tests whether the observed mean across flies exceeds this null (one-sided).


nDraw_pop = 100000; 

% Slope
has_sn_matched = cellfun(@(x) ~isempty(x), S.slope_dist_sn_matched_pooled_all);
has_actual_sn  = isfinite(S.slope_actual_all);
valid_pop_sn   = has_sn_matched & has_actual_sn;

if sum(valid_pop_sn) < numel(S.fly_IDs_all)
    warning('%d of %d flies excluded from shot noise test', ...
        numel(S.fly_IDs_all) - sum(valid_pop_sn), numel(S.fly_IDs_all));
end

fly_nulls_sn_matched = S.slope_dist_sn_matched_pooled_all(valid_pop_sn);
nFlies_sn_matched    = numel(fly_nulls_sn_matched);
obs_mean_slope_sn    = mean(S.slope_actual_all(valid_pop_sn), 'omitnan');

pop_null_sn_matched = nan(nDraw_pop, 1);

for b = 1:nDraw_pop
    sampled_sn = cellfun(@(v) v(randi(numel(v))), fly_nulls_sn_matched);
    pop_null_sn_matched(b) = mean(sampled_sn);
end

threshold_95_sn   = prctile(pop_null_sn_matched, 95); % cutoff, NOT a confidence interval
null_median_sn    = median(pop_null_sn_matched);
null_mean_sn      = mean(pop_null_sn_matched);

n_exc_pop_sn_slope = sum(pop_null_sn_matched >= obs_mean_slope_sn);
p_pop_sn_matched   = (1 + n_exc_pop_sn_slope) / (1 + nDraw_pop);

fprintf('\n--- Population matched-sample shot noise test: slope ---\n');
fprintf('  N flies:                        %d\n',   nFlies_sn_matched);
fprintf('  Observed mean actual slope:     %.4f\n', obs_mean_slope_sn);
fprintf('  Null median:                    %.4f\n', null_median_sn);
fprintf('  Null mean:                      %.4f\n', null_mean_sn);
fprintf('  95th percentile threshold:      %.4f\n', threshold_95_sn);
fprintf('  p (one-sided, Monte Carlo):     %.2g\n', p_pop_sn_matched);
fprintf('  %d/%d null >= obs\n', n_exc_pop_sn_slope, nDraw_pop);

S.pop.sn_matched_slope.obs_mean       = obs_mean_slope_sn;
S.pop.sn_matched_slope.null_dist      = pop_null_sn_matched;
S.pop.sn_matched_slope.p_pop          = p_pop_sn_matched;
S.pop.sn_matched_slope.threshold_95   = threshold_95_sn;
S.pop.sn_matched_slope.null_median    = null_median_sn;
S.pop.sn_matched_slope.null_mean      = null_mean_sn;
S.pop.sn_matched_slope.nFlies         = nFlies_sn_matched;
S.pop.sn_matched_slope.n_exceed       = n_exc_pop_sn_slope;
S.pop.sn_matched_slope.nDraw          = nDraw_pop;

% Rho

has_sn_rho     = cellfun(@(x) ~isempty(x), S.rho_dist_sn_matched_pooled_all);
valid_pop_sn_rho = has_sn_rho & isfinite(S.rho_actual_all);

if sum(valid_pop_sn_rho) < numel(S.fly_IDs_all)
    warning('%d of %d flies excluded from population shot noise rho test.', ...
        numel(S.fly_IDs_all) - sum(valid_pop_sn_rho), numel(S.fly_IDs_all));
end

fly_nulls_sn_rho = S.rho_dist_sn_matched_pooled_all(valid_pop_sn_rho);
nFlies_sn_rho = numel(fly_nulls_sn_rho);
obs_mean_rho_sn  = mean(S.rho_actual_all(valid_pop_sn_rho), 'omitnan');

pop_null_sn_rho = nan(nDraw_pop, 1);
for b = 1:nDraw_pop
    pop_null_sn_rho(b) = mean(cellfun(@(v) v(randi(numel(v))), fly_nulls_sn_rho));
end

threshold_95_sn_rho = prctile(pop_null_sn_rho, 95);
null_median_sn_rho = median(pop_null_sn_rho);
null_mean_sn_rho = mean(pop_null_sn_rho);

n_exc_pop_sn_rho = sum(pop_null_sn_rho >= obs_mean_rho_sn);
p_pop_sn_rho     = (1 + n_exc_pop_sn_rho) / (1 + nDraw_pop);

fprintf('\n--- Population matched-sample shot noise test: rho ---\n');
fprintf('  N flies:                        %d\n',  nFlies_sn_rho);
fprintf('  Observed mean actual rho:       %.4f\n', obs_mean_rho_sn);
fprintf('  Null median:                    %.4f\n', null_median_sn_rho);
fprintf('  Null mean:                      %.4f\n', null_mean_sn_rho);
fprintf('  95th percentile threshold:      %.4f\n', threshold_95_sn_rho);
fprintf('  p (one-sided, Monte Carlo):     %.2g\n', p_pop_sn_rho);
fprintf('  %d/%d null >= obs\n', n_exc_pop_sn_rho, nDraw_pop);


S.pop.sn_matched_rho.obs_mean     = obs_mean_rho_sn;
S.pop.sn_matched_rho.null_dist    = pop_null_sn_rho;
S.pop.sn_matched_rho.p_pop        = p_pop_sn_rho;
S.pop.sn_matched_rho.threshold_95 = threshold_95_sn_rho;
S.pop.sn_matched_rho.null_median  = null_median_sn_rho;
S.pop.sn_matched_rho.null_mean    = null_mean_sn_rho;
S.pop.sn_matched_rho.nFlies       = nFlies_sn_rho;
S.pop.sn_matched_rho.n_exceed     = n_exc_pop_sn_rho;
S.pop.sn_matched_rho.nDraw        = nDraw_pop;



% ---------------- SAVE EVERYTHING ----------------

save(summaryFile, '-struct', 'S');

fprintf('Master summary updated with population statistics:\n%s\n', summaryFile);



%% ------------------------------------------------------------------------
%                            4. PLOTTING & SAVING
%
%                Creates one 'Dashboard' window with tabs
% -------------------------------------------------------------------------

fprintf('\n--- Generating and Saving Figures ---\n');


% --- 17) SETUP DASHBOARD WINDOW --- 
% Check if dashboard already exists to prevent duplicates

dashboardName = 'MBON Analysis Dashboard';
hDash = findall(0, 'Type', 'figure', 'Name', dashboardName);

if isempty(hDash)
    hDash = figure('Name', dashboardName, 'Color', 'w', 'NumberTitle', 'off', ...
                   'Position', [50, 50, 1200, 800]); 
    tabGroup = uitabgroup(hDash); % Create the tab group container
else
    figure(hDash); % Bring to front
    tabGroup = findobj(hDash, 'Type', 'uitabgroup'); % Find existing tabgroup to append to
    if isempty(tabGroup)
        tabGroup = uitabgroup(hDash); % Create if missing
    end
end


% 18) --- Dashboard Calls: Save, then dock ---
% This helper function handles the "Save first, then Tab (Display)" logic,
% which preserves resolution in .svg 

% Calls plotting helper functions
% Uses 'Anonymous Functions' (@()...) to bundle variables with functions
% before passsing to the helper. This is needed because helper functions
% cannot access main workspace variables

% --- Per-fly DF/F traces ---
add_plot_to_dashboard(tabGroup, [flyID ' - Mean Traces'], ...
    fullfile(figsFolder, [flyID '_MeanTraces.svg']), @() plot_mean_traces(...
    all_traces_ROI1, all_traces_ROI2, stimFrames, flyID, fs));

add_plot_to_dashboard(tabGroup, [flyID ' - Amplitude vs Frequency'], ...
    fullfile(figsFolder, [flyID '_AmplitudeFrequency.svg']), ...
    @() plot_amplitude_frequency(all_amplitudes_ROI1, all_amplitudes_ROI2, flyID, frequency), ...
    [800*2/3, 600]);

% --- Per-fly TLS regression  ---
add_hybrid_plot_to_dashboard(tabGroup, [flyID ' - Shot MBONap TLS QC'], ...
    fullfile(figsFolder, [flyID '_Shot_ap_TLS_QC_data.png']), ...
    fullfile(figsFolder, [flyID '_Shot_ap_TLS_QC_frame.pdf']), ...
    @() plot_tls_qc(A1_long, B1_long, v_tls_s1, proj_sign_s1, r_orth_s1_scaled, flyID, 'Shot MBONap', false, ...
        'ROI1a amplitude', 'ROI1b amplitude', shared_xlim, shared_ylim, 'data_only'), ...
    @() plot_tls_qc(A1_long, B1_long, v_tls_s1, proj_sign_s1, r_orth_s1_scaled, flyID, 'Shot MBONap', false, ...
        'ROI1a amplitude', 'ROI1b amplitude', shared_xlim, shared_ylim, 'frame_only'), ...
    @() plot_tls_qc(A1_long, B1_long, v_tls_s1, proj_sign_s1, r_orth_s1_scaled, flyID, 'Shot MBONap', false, ...
        'ROI1a amplitude', 'ROI1b amplitude', shared_xlim, shared_ylim, 'full'));

add_hybrid_plot_to_dashboard(tabGroup, [flyID ' - Shot MBONav TLS QC'], ...
    fullfile(figsFolder, [flyID '_Shot_av_TLS_QC_data.png']), ...
    fullfile(figsFolder, [flyID '_Shot_av_TLS_QC_frame.pdf']), ...
    @() plot_tls_qc(A2_long, B2_long, v_tls_s2, proj_sign_s2, r_orth_s2_scaled, flyID, 'Shot MBONav', false, ...
        'ROI2a amplitude', 'ROI2b amplitude', shared_xlim, shared_ylim, 'data_only'), ...
    @() plot_tls_qc(A2_long, B2_long, v_tls_s2, proj_sign_s2, r_orth_s2_scaled, flyID, 'Shot MBONav', false, ...
        'ROI2a amplitude', 'ROI2b amplitude', shared_xlim, shared_ylim, 'frame_only'), ...
    @() plot_tls_qc(A2_long, B2_long, v_tls_s2, proj_sign_s2, r_orth_s2_scaled, flyID, 'Shot MBONav', false, ...
        'ROI2a amplitude', 'ROI2b amplitude', shared_xlim, shared_ylim, 'full'));

add_plot_to_dashboard(tabGroup, [flyID ' - Actual MBON TLS QC'], ...
    fullfile(figsFolder, [flyID '_Actual_TLS_QC.svg']), ...
    @() plot_tls_qc(amps1, amps2, v_tls_actual, proj_sign_actual, r_orth_actual, flyID, 'Actual MBON', true, ...
    [], [], shared_xlim, shared_ylim));

% --- Per-fly OLS ---
add_plot_to_dashboard(tabGroup, [flyID ' - Noise Scaling'], ...
    fullfile(figsFolder, [flyID '_NoiseScaling.svg']), ...
    @() plot_ols_noise_scaling(...
        proj_sign_s1, r_orth_s1_scaled, ...
        proj_sign_s2, r_orth_s2_scaled, ...
        proj_sign_actual, r_orth_actual, ...
        flyID));

% --- Plot Permutation Null ----
add_plot_to_dashboard(tabGroup, [flyID ' - Permutation Null'], ...
    fullfile(figsFolder, [flyID '_PermNull.svg']), ...
    @() plot_permutation_null(null_rho, null_slope, ...
        rho_actual, slope_act, p_perm_rho_actual, ...
        p_perm_slope_actual, flyID));

% --- Per-fly: matched shot noise slope distributions vs actual slope ---
add_plot_to_dashboard(tabGroup, [flyID ' - Matched SN Control'], ...
    fullfile(figsFolder, [flyID '_MatchedSN_Control.svg']), ...
    @() plot_sn_matched_control( ...
        slope_sn_matched_pooled, rho_sn_matched_pooled, ...
        slope_act, rho_actual, ...
        p_sn_matched_slope, p_sn_matched_rho, flyID));

% --- Population summary - timestamped so each run's snapshot is preserved ---
popTimestamp     = char(datetime('now', 'Format', 'yyyy-MM-dd_HHmmss'));
populationFigPath = fullfile(figsFolder, ['Population_Summary_Deming_WithinVar_' popTimestamp '.svg']);

add_plot_to_dashboard(tabGroup, 'Population Summary', populationFigPath, ...
    @() plot_population_summary(S));

add_plot_to_dashboard(tabGroup, 'Population - Slope Null Test', ...
    fullfile(figsFolder, ['Population_SlopeNull_' popTimestamp '.svg']), ...
    @() plot_slope_null_population(S));

add_plot_to_dashboard(tabGroup, 'Population - Rho Null Test', ...
    fullfile(figsFolder, ['Population_RhoNull_' popTimestamp '.svg']), ...
    @() plot_rho_null_population(S));

add_plot_to_dashboard(tabGroup, 'Population - All Flies Signal Variability', ...
    fullfile(figsFolder, ['Population_AllFlies_SignalVariability_' popTimestamp '.svg']), ...
    @() plot_population_signal_variability(S), [400 600]);

add_plot_to_dashboard(tabGroup, 'Population - Matched SN Null Test', ...
    fullfile(figsFolder, ['Population_MatchedSN_' popTimestamp '.svg']), ...
    @() plot_sn_matched_population(S));

% --- Combined violin plots: population + per-fly null distributions ---
add_plot_to_dashboard(tabGroup, 'Population - Permutation Null Violins', ...
    fullfile(figsFolder, ['Population_PermNull_Violins_' popTimestamp '.svg']), ...
    @() plot_combined_null_violins_perm(S));

add_plot_to_dashboard(tabGroup, 'Population - Matched SN Violins', ...
    fullfile(figsFolder, ['Population_MatchedSN_Violins_' popTimestamp '.svg']), ...
    @() plot_combined_null_violins_sn(S));

% --- Save dashboard as .fig - timestamped to preserve each run's snapshot ---
dashboardFile = fullfile(figsFolder, ['MBON_Dashboard_Deming_WithinVar_' flyID '_' popTimestamp '.fig']);
savefig(hDash, dashboardFile);
fprintf('Dashboard saved to: %s\n', dashboardFile);
fprintf('All figures saved to: %s\n', figsFolder);

diary off;


%% ----------------------------------------------------------------------------
%                        HELPER FUNCTIONS - PLOTTING
%
%                  Define *how* to plot - do not open windows
% -----------------------------------------------------------------------------
% NB: number of input arguments needs to match those in call (unless
% varargin), BUT for flexibility need to rename local input arguments and not use global
% workspace variables


% --- Handles the "Save Invisible -> Show in Tab" logic ---
function add_plot_to_dashboard(tg, tabTitle, savePath, plotFunc, figSize)
    if nargin < 5 || isempty(figSize)
        figSize = [800 600]; % default, for two panel plots
    end
    % 1. Create invisible figure for saving
    hTemp = figure('Visible', 'off', 'Color', 'w', 'Position', [0 0 figSize]);

   % 2. Set painters renderer for clean vector output and editability in Illustrator 
     % Rendering with opengl often rasterises/outlines objects during
        % export and converts text into individual vector paths, rather than
        % keeping it as editable text
     %set(hTemp, 'Renderer', 'painters'); % 'renderer' will be removed in future release, Matlab switching from openGL to WebGL -> Add -vector flag to print 

    % 3. Run the plotting function
    plotFunc();

    % 4. Save as SVG (Illustrator editable)
     print(hTemp, savePath, '-dsvg', '-vector');
    % exportgraphics(hTemp, savePath, 'ContentType', 'vector', 'BackgroundColor', 'white'); % NB: does not work with special/greek characters, rendered by tex interpreter, 
                                                                                            % all text ends up as outlines paths
    % 5. Create Tab
    t = uitab(tg, 'Title', tabTitle);

    % 6. Copy content (explicitly add axes and legends, instead of
    % previous, unstable 'copyobj(allchild(hTemp),t)'
    objs = findobj(hTemp, 'Type', 'axes', '-or', 'Type', 'legend'); % expand if other items such as colorbar needed
    copyobj(objs, t);

    close(hTemp);
end


% --- Handles hybrid raster (data) + vector (axes/labels) export for dense scatter plots ---
function add_hybrid_plot_to_dashboard(tg, tabTitle, savePathRaster, savePathVector, plotFuncRaster, plotFuncVector, plotFuncFull, figSizeInches)
    if nargin < 8 || isempty(figSizeInches)
        figSizeInches = [8 6];   % width, height in inches
    end
    
    % 1. Raster: data only
    hRaster = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
                      'Position', [0 0 figSizeInches]);
    plotFuncRaster();
    print(hRaster, savePathRaster, '-dpng', '-r600');
    close(hRaster);

    % 2. Vector: axes frame, ticks, labels only
    hVector = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
                      'Position', [0 0 figSizeInches]);
    set(hVector, 'PaperUnits', 'inches', ...
                 'PaperPositionMode', 'manual', ...
                 'PaperPosition', [0 0 figSizeInches], ...
                 'PaperSize', figSizeInches);
    plotFuncVector();
    print(hVector, savePathVector, '-dpdf', '-vector');
    close(hVector);

    % 3. Full: for dashboard tab QC viewing only
    hFull = figure('Visible', 'off', 'Color', 'w', 'Units', 'inches', ...
                    'Position', [0 0 figSizeInches]);
    plotFuncFull();
    t = uitab(tg, 'Title', tabTitle);
    objs = findobj(hFull, 'Type', 'axes', '-or', 'Type', 'legend');
    copyobj(objs, t);
    close(hFull);
end


% --- Define CVD-friendly colormap - Okabe Ito ---

function cmap = okabe_ito_colors()
    cmap.black      = [0 0 0];
    cmap.orange     = [0.90 0.60 0.00];
    cmap.skyblue    = [0.35 0.70 0.90];
    cmap.green      = [0.00 0.60 0.50];
    cmap.yellow     = [0.95 0.90 0.25];
    cmap.blue       = [0.00 0.45 0.70];
    cmap.vermillion = [0.80 0.40 0.00];
    cmap.purple     = [0.80 0.60 0.70];
end


% --- Display all traces w. mean, plus per-frequency amplitude panel ---
function plot_mean_traces(traces1, traces2, stim, id, fs)
% traces1, traces2  : [nFrames x nTrials] data
% stim              : stimulus onset frames
% id                : identifier string
% fs                : sampling rate in Hz (frames per second)

        c = okabe_ito_colors();
        nFrames = size(traces1, 1);
        t = (0:nFrames-1) / fs;   % time vector in seconds
        stim_sec = stim / fs;     % convert stim onset frames to seconds 

% --- Panel 1: ROI1 traces (wide, left) ---
        subplot(2,1,1);
        set(gca, 'TickDir', 'out');
        hold on;
        mean_trace_roi1 = mean(traces1, 2, 'omitnan');
        plot(t, traces1, 'Color', [0.7 0.7 0.7, 0.5], 'LineWidth', 0.5);
        plot(t, mean_trace_roi1, 'k', 'LineWidth', 2.0);
        for k = 1:numel(stim_sec)
            xline(stim_sec(k), '--', 'Color', c.vermillion, 'LineWidth', 1);
        end
        title(['All Traces ROI 1: ' id], 'Interpreter', 'none');
        ylabel('dF/F_0');
        xlim([t(1) - range(t)*0.02, t(end) + range(t)*0.02]);
        box off;

% --- Panel 2: ROI2 traces (wide, left) ---
        subplot(2, 1, 2);
        set(gca, 'TickDir', 'out');
        hold on;
        mean_trace_roi2 = mean(traces2, 2, 'omitnan');
        plot(t, traces2, 'Color', [0.7 0.7 0.7, 0.5], 'LineWidth', 0.5);
        plot(t, mean_trace_roi2, 'k', 'LineWidth', 2.0);
        for k = 1:numel(stim_sec)
            xline(stim_sec(k), '--', 'Color', c.vermillion, 'LineWidth', 1);
        end
        title(['All Traces ROI 2: ' id], 'Interpreter', 'none');
        xlabel('Time (s)');
        ylabel('dF/F_0');
        xlim([t(1) - range(t)*0.02, t(end) + range(t)*0.02]);
        box off;
end


% --- Display amplitude vs frequency (ROI1 + ROI2) ---
function plot_amplitude_frequency(amps1, amps2, id, freqs)
% amps1, amps2 : [nStim x nTrials] peak amplitude per stimulus/trial
% id           : flyID
% freqs        : [nStim x 1] stimulus frequencies, e.g. [2 4 10 20 40]

        if nargin < 4 || isempty(freqs)
            freqs = [2;4;10;20;40];
        end
        freqs = freqs(:);

% --- Panel 1: ROI1 amplitude vs frequency ---
        subplot(2, 1, 1);
        set(gca, 'TickDir', 'out');
        hold on;
        yline(0, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        mean_amp_roi1 = mean(amps1, 2, 'omitnan');
        plot(freqs, amps1, 'Color', [0.7 0.7 0.7, 0.5], 'LineWidth', 0.3);
        scatter(repmat(freqs, 1, size(amps1,2)), amps1, 25, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.7);
        scatter(freqs, mean_amp_roi1, 40, 'k', 'filled');
        set(gca, 'XScale', 'log', 'XTick', freqs, 'XTickLabel', freqs);
        xlim([freqs(1) / 1.3, freqs(end) * 1.3]);
        xlabel('Frequency (Hz, log scale)');
        ylabel('Peak dF/F_0');
        title('ROI 1 amplitude');
        box off;

% --- Panel 2: ROI2 amplitude vs frequency ---
        subplot(2, 1, 2);
        set(gca, 'TickDir', 'out');
        hold on;
        yline(0, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        mean_amp_roi2 = mean(amps2, 2, 'omitnan');
        plot(freqs, amps2, 'Color', [0.7 0.7 0.7, 0.5], 'LineWidth', 0.3);
        scatter(repmat(freqs, 1, size(amps2,2)), amps2, 25, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.7);
        scatter(freqs, mean_amp_roi2, 40, 'k', 'filled');
        set(gca, 'XScale', 'log', 'XTick', freqs, 'XTickLabel', freqs);
        xlim([freqs(1) / 1.3, freqs(end) * 1.3]);
        xlabel('Frequency (Hz, log scale)');
        ylabel('Peak dF/F_0');
        title('ROI 2 amplitude');
        box off;
end


 
% --- Generic helper function for TLS QC (use w. diff inputs) - two
% panels: left scatter with TLS axis, right r_orth vs tls_projection
function plot_tls_qc(A, B, v_tls, proj_sign, r_orth, flyID, datasetName, draw_lines, xlab, ylab, shared_xlim, shared_ylim, mode)
    % TLS QC plot per fly / dataset
    % A, B           - response vectors
    % v_tls          - TLS unit vector
    % proj_sign      - signed projection along TLS axis
    % r_orth         - orthogonal residual magnitude
    % flyID          - fly identifier string
    % datasetName    - dataset name string
    % draw_lines     - boolean, whether to draw perpendicular lines
    % xlab, ylab     - optional panel-1 axis labels
    % shared_x/ylim  - optional shared axis limits for panel 2
    % mode           - 'full' (default), 'data_only', or 'frame_only'
    %                  'data_only'  : plot scatter points only, no axes/labels/ticks
    %                  'frame_only' : plot axes, ticks, labels, TLS line, title only, no scatter clouds
    %                  Used to split dense scatter data (rasterized) from
    %                  editable vector frame/text (kept as vector) for
    %                  Illustrator-friendly export of high point-count figures

    if nargin < 8
        draw_lines = false;
    end
    if nargin < 9 || isempty(xlab)
        xlab = 'MBONap amplitude';
    end
    if nargin < 10 || isempty(ylab)
        ylab = 'MBONav amplitude';
    end
    if nargin < 11
        shared_xlim = [];
    end
    if nargin < 12
        shared_ylim = [];
    end
    if nargin < 13 || isempty(mode)
        mode = 'full';
    end

    plot_data  = ~strcmp(mode, 'frame_only');
    plot_frame = ~strcmp(mode, 'data_only');

    % Set colormap
    c = okabe_ito_colors();

    % Input
    X = [A(:), B(:)];

    % --- Panel 1: ROI1 vs ROI2 with TLS axis ---
    subplot(1,2,1);
    hold on;

    if plot_data
        if draw_lines
            scatter(A, B, 25, c.black, 'filled', 'MarkerEdgeColor', 'none');
        else
            %scatter(A, B, 25, 'MarkerEdgeColor', c.black, 'MarkerEdgeAlpha', 0.1, 'LineWidth', 0.5);
            scatter(A, B, 25, c.black, 'filled', 'MarkerFaceAlpha', 0.05);
        end
    end

    % --- STEP 1: Define the anchor (Mathematical Centroid) ---
    mu = mean(X, 1);

    % --- STEP 2: Determine draw extents based on data range ---
    centered_X = X - mu;
    projections = centered_X * v_tls(:);

    t_min = min(projections);
    t_max = max(projections);
    buffer = (t_max - t_min) * 0.1;
    t_start = t_min - buffer;
    t_end   = t_max + buffer;

    pt_start = mu + t_start * v_tls(:)';
    pt_end   = mu + t_end   * v_tls(:)';

    if plot_frame
        plot([pt_start(1) pt_end(1)], [pt_start(2) pt_end(2)], 'Color', c.vermillion, 'LineWidth', 2);
        xlabel(xlab);
        ylabel(ylab);
        title([datasetName ' - TLS fit']);
    end

    % --- STEP 3: Equal axes (Square Aspect Ratio) ---
    axis equal;
    global_min = min(min(X));
    global_max = max(max(X));
    range_buffer = (global_max - global_min) * 0.05;
    lims = [global_min - range_buffer, global_max + range_buffer];
    xlim(lims);
    ylim(lims);

    if plot_frame
        tick_start  = floor(lims(1) / 0.5) * 0.5;
        tick_end    = ceil(lims(2) / 0.5) * 0.5;
        ticks_05    = tick_start:0.5:tick_end;
        set(gca, 'XTick', ticks_05, 'YTick', ticks_05);
        set(gca, 'TickDir', 'out');
        box off;
    else
        set(gca, 'Visible', 'off');
    end

    % --- Draw perpendicular lines (Residuals) ---
    if draw_lines && plot_data
        pts_on_axis = mu + projections * v_tls(:)';
        for i = 1:size(X,1)
             line([X(i,1), pts_on_axis(i,1)], [X(i,2), pts_on_axis(i,2)], ...
                'Color', [c.black, 0.5], 'LineWidth', 0.5);
        end
    end

    % --- Panel 2: Projection vs orthogonal residual magnitude ---
    subplot(1,2,2);
    hold on;

    if plot_data
        if draw_lines
            scatter(proj_sign, r_orth, 25, c.black, 'filled', 'MarkerFaceAlpha',1);
        else
            scatter(proj_sign, r_orth, 25, c.black, 'filled', 'MarkerFaceAlpha',0.05);
        end
    end

    if plot_frame
        xlabel('Projection onto TLS (shifted)');
        ylabel('Orthogonal distance');
    end

    axis square;
    if ~isempty(shared_xlim)
        xlim(shared_xlim);
    end
    if ~isempty(shared_ylim)
        ylim(shared_ylim);
    end

    % --- Shot noise slope fit (only for shot noise, not actual data) ---
    if ~draw_lines && plot_data
        valid = isfinite(proj_sign) & isfinite(r_orth);
        mdl_sn = fitlm(proj_sign(valid), r_orth(valid));
        x_sn_fit = linspace(min(proj_sign(valid)), max(proj_sign(valid)), 100)';
        y_sn_fit = mdl_sn.Coefficients.Estimate(1) + mdl_sn.Coefficients.Estimate(2) * x_sn_fit;
        plot(x_sn_fit, y_sn_fit, 'Color', c.skyblue, 'LineWidth', 2);
    end

    if plot_frame
        set(gca, 'TickDir', 'out');
        box off;

        % Spearman correlation
        [rho, pval] = corr(proj_sign, r_orth, 'Type', 'Spearman', 'Rows','complete');
        title(sprintf('%s\nSpearman rho=%.3f, p=%.3f', datasetName, rho, pval));
    else
        set(gca, 'Visible', 'off');
    end

    if plot_frame
        sgtitle(sprintf('Fly %s: %s TLS QC', flyID, datasetName));
    end
end

% ---
function plot_ols_noise_scaling(proj1, rorth1, proj2, rorth2, projAct, rorthAct, flyID)

    % Per-fly Noise Scaling: Linear OLS 
    datasets = {'Shot MBONap','Shot MBONav','Actual'};
    proj_all = {proj1, proj2, projAct};
    rorth_all = {rorth1, rorth2, rorthAct};
    colors = {[0.8 0.2 0.2],[0.2 0.2 0.8],[0 0 0]};
    
    for i = 1:3
        x = proj_all{i};
        y = rorth_all{i};
        
        subplot(1,3,i); 
        set(gca, 'TickDir', 'out');
        hold on;

        valid = isfinite(x) & isfinite(y);
        x_clean = x(valid);
        y_clean = y(valid);

        scatter(x_clean, y_clean, 25, colors{i}, 'filled', 'MarkerFaceAlpha',0.5);
        mdl = fitlm(x_clean, y_clean);

        % For red line 
        x_fit = linspace(min(x_clean), max(x_clean), 100)';
        y_fit = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2) * x_fit;
        plot(x_fit, y_fit, 'r-', 'LineWidth', 2);

        xlabel('Projection (latent drive)');
        ylabel('Orthogonal distance (noise)');
        title([datasets{i} ' - Linear slope: ' num2str(mdl.Coefficients.Estimate(2),'%.2f')]);
        axis square;
        box off;
    end
    
    sgtitle(['Fly ' flyID ': Linear Noise Scaling']);
end


    % ---
function plot_population_summary(S)
    labels = {'Shot MBONap', 'Shot MBONav', 'Actual'};
    colors = {[0.8 0.2 0.2], [0.2 0.2 0.8], [0 0 0]};

    % Set colormap
    c = okabe_ito_colors();
    
    % Use a layout factor (0.66 = 2/3 width)
    width_scale = 0.66; 

    % --- SPEARMAN RHO (Top Subplot) ---
    ax1 = make_narrow_subplot(2, 1, 1, width_scale); 
    hold on;
    groups_rho = {S.rho_shot1_all, S.rho_shot2_all, S.rho_actual_all};
    ymax = -inf; 
    ymin = inf;
    
    for g = 1:3
        r = groups_rho{g};
        r = r(isfinite(r));
        if isempty(r)
            continue; 
        end
        scatter(g*ones(size(r)), r, 35, colors{g}, 'filled', ...
            'MarkerFaceAlpha',0.7, 'Jitter','on', 'JitterAmount',0.1);
        plot([g-0.25 g+0.25], [mean(r) mean(r)], 'Color', c.vermillion, 'LineWidth', 2); % plot mean instead of median here, empirical p_value (permutation null) is now main reported result rather than wilcoxon against 0
        ymax = max(ymax, max(r));
        ymin = min(ymin, min(r));
    end
    yline(0, '--', 'Color', [0.6 0.6 0.6]);
    set(gca, 'XTick', 1:3, 'XTickLabel', labels);
    ylabel('Spearman rho');
    title('Signal dependent noise');
    
    % Build Text for Right Side (Wilcoxon vs 0)
    txt_rho = '\bf Wilcoxon vs 0: \rm'; % Header
    group_names = {'Shot1', 'Shot2', 'Actual'};
    if isfield(S.pop, 'rho_wilcoxon')
        for d = 1:min(3, numel(S.pop.rho_wilcoxon))
            p = S.pop.rho_wilcoxon(d).p_value;
            txt_rho = [txt_rho, sprintf('\n%s: %.3g (%s)', group_names{d}, p, pstars(p))];
        end
    end
    
    
    % Add paired comparisons
    if isfield(S.pop, 'rho_wilcoxon_paired') && numel(S.pop.rho_wilcoxon_paired) >= 2 % 2 refers to comparisons (actual vs shot1, actuak vs shot2)
        txt_rho = [txt_rho, sprintf('\n\n\\bf Paired Wilcoxon (vs Act): \\rm')];
        p1 = S.pop.rho_wilcoxon_paired(1).p_value;
        p2 = S.pop.rho_wilcoxon_paired(2).p_value;
        txt_rho = [txt_rho, sprintf('\nShot1: %.3g (%s)', p1, pstars(p1))];
        txt_rho = [txt_rho, sprintf('\nShot2: %.3g (%s)', p2, pstars(p2))];
    end
    
    % Paired permutation 
    if isfield(S.pop, 'rho_paired_perm') 
        txt_rho = [txt_rho, sprintf('\n\n\\bf Paired permutation: \\rm')];
        txt_rho = [txt_rho, sprintf('\nShot1: %.3g', S.pop.rho_paired_perm.p_shot1)];
        txt_rho = [txt_rho, sprintf('\nShot2: %.3g', S.pop.rho_paired_perm.p_shot2)];
    end
    
    text(1.05, 0.5, txt_rho, 'Units', 'normalized', ...
        'VerticalAlignment', 'middle', 'FontSize', 9, 'Color', 'k');

    beautify_axis();
    box off;


    % --- OLS SLOPES (Bottom Subplot) ---
    ax2 = make_narrow_subplot(2, 1, 2, width_scale); 
    hold on;
    groups_slope = {S.slope_shot1_all, S.slope_shot2_all, S.slope_actual_all};
    ymax = -inf; ymin = inf;
    
    for g = 1:3
        r = groups_slope{g};
        r = r(isfinite(r));
        if isempty(r) 
            continue; 
        end
        scatter(g*ones(size(r)), r, 35, colors{g}, 'filled', ...
            'MarkerFaceAlpha',0.7, 'Jitter','on', 'JitterAmount',0.1);
        plot([g-0.25 g+0.25], [mean(r) mean(r)], 'Color', c.vermillion, 'LineWidth', 2); % mean instead of median, to match permutation null control
        ymax = max(ymax, max(r));
        ymin = min(ymin, min(r));
    end
    yline(0, '--', 'Color', [0.6 0.6 0.6]);
    set(gca, 'XTick', 1:3, 'XTickLabel', labels);
    ylabel('OLS slope');
    title('Linear scaling across flies');
    
    % Build Combined Text for Right Side
    % Paired Wilcoxon
    txt_slope = '\bf Paired Wilcoxon (vs Act): \rm';
    if isfield(S.pop, 'slope_wilcoxon') && numel(S.pop.slope_wilcoxon) >= 2
        p1 = S.pop.slope_wilcoxon(1).p_value;
        p2 = S.pop.slope_wilcoxon(2).p_value;
        txt_slope = [txt_slope, sprintf('\nShot1: %.3g (%s)', p1, pstars(p1))];
        txt_slope = [txt_slope, sprintf('\nShot2: %.3g (%s)', p2, pstars(p2))];
    end
    
    % Paired Permutation (append if exists)
    if isfield(S.pop, 'slope_paired_perm') 
        txt_slope = [txt_slope, sprintf('\n\n\\bf Paired permutation: \\rm')];
        txt_slope = [txt_slope, sprintf('\nShot1: %.3g', S.pop.slope_paired_perm.p_shot1)];
        txt_slope = [txt_slope, sprintf('\nShot2: %.3g', S.pop.slope_paired_perm.p_shot2)];
    end
    
    text(1.05, 0.5, txt_slope, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 9, 'Color', [0.2 0.2 0.2]); % Dark grey for readability

    beautify_axis();
    box off;

    sgtitle('Summary Across Flies');
end


% --- 
function plot_permutation_null(null_rho, null_slope, ...
        rho_actual, slope_act, ...
        p_perm_rho, p_perm_slope, flyID)
% Per-fly permutation null plot: two panels
%   Left   : Spearman rho distribution
%   Right  : OLS slope distribution (raw, uncorrected)
% Observed actual value shown as red line; null median as dashed black line.

    panel_data = { ...
        null_rho,   rho_actual, p_perm_rho,   'Spearman rho'; ...
        null_slope, slope_act,  p_perm_slope,  'OLS slope (raw)'};
 
    for i = 1:2
        subplot(1, 2, i);
        hold on;
 
        dist     = panel_data{i, 1};
        obs_val  = panel_data{i, 2};
        p_val    = panel_data{i, 3};
        stat_lbl = panel_data{i, 4};
 
        histogram(dist, 50, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', ...
            'DisplayName', 'Null dist');
 
        yl = ylim;
 
        plot([obs_val obs_val], yl, 'r-', 'LineWidth', 2.5, ...
            'DisplayName', sprintf('Observed (%.3f)', obs_val));
 
        med_val = median(dist, 'omitnan');
        plot([med_val med_val], yl, 'k--', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Median (%.3f)', med_val));
 
        xlabel(stat_lbl);
        ylabel('Count');
        title(sprintf('Null distribution: %s\np_{perm} = %.4f', stat_lbl, p_val), ...
            'Interpreter', 'tex');
        legend('Location', 'northwest', 'FontSize', 8);
        box off;
        set(gca, 'TickDir', 'out');
        box off;
    end
 
    sgtitle(sprintf('Fly %s: Permutation null (raw slopes)', flyID), ...
        'Interpreter', 'none');
end


function plot_rho_null_population(S)
% Population-level permutation null rho test
% Monte Carlo null distribution with observed mean rho
% NB: null is centred ABOVE zero due to TLS geometric bias

    null_dist     = S.pop.rho_perm.null_dist;
    obs_mean      = S.pop.rho_perm.obs_mean;
    null_mean     = S.pop.rho_perm.null_mean;
    p_val         = S.pop.rho_perm.p_pop_rho;
    threshold_95  = S.pop.rho_perm.threshold_95;
    nFlies        = S.pop.rho_perm.nFlies;


    colors.null     = [0.7 0.7 0.7];
    colors.observed = [0.8 0.2 0.2];

    hold on;

    histogram(null_dist, 60, ...
        'FaceColor', colors.null, 'EdgeColor', 'none', ...
        'DisplayName', 'Permutation null (TLS geometry)');

    yl = ylim;

    % 95th percentile shading
    x_shade = [threshold_95, max(null_dist), ...
                max(null_dist), threshold_95];
    y_shade = [0, 0, yl(2), yl(2)];
    fill(x_shade, y_shade, [0.9 0.7 0.7], ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
        'DisplayName', '95th percentile');

    % Null mean (the null is a null of per-fly means, so compare means for consistency)
    plot([null_mean null_mean], yl, 'k--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Null mean (%.3f)', null_mean));

    % Observed mean rho
    plot([obs_mean obs_mean], yl, '-', ...
        'Color', colors.observed, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Observed mean rho (%.3f)', obs_mean));

    % Zero reference
    plot([0 0], yl, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
        'HandleVisibility', 'off');

    xlabel('Mean Spearman rho');
    ylabel('Count');
    title(sprintf(['Population permutation null (rho)\n' ...
        'Null mean = %.3f (above zero: TLS geometric bias)\n' ...
        'p = %.4f (n = %d flies)'], ...
        null_mean, p_val, nFlies));
    legend('Location', 'northwest', 'FontSize', 8);
    box off;
    set(gca, 'TickDir', 'out');
end

% ---
function plot_slope_null_population(S)
% Population-level permutation null slope test
% Monte Carlo null distribution with observed mean slope

    null_dist     = S.pop.slope_perm.null_dist;
    obs_mean      = S.pop.slope_perm.obs_mean;
    null_mean     = S.pop.slope_perm.null_mean;
    p_val         = S.pop.slope_perm.p_pop_slope;
    threshold_95  = S.pop.slope_perm.threshold_95;
    nFlies        = S.pop.slope_perm.nFlies;
   
    colors.null     = [0.7 0.7 0.7];
    colors.observed = [0.8 0.2 0.2];
 
 
    hold on;
 
    histogram(null_dist, 60, ...
        'FaceColor', colors.null, 'EdgeColor', 'none', ...
        'DisplayName', 'Permutation null');
 
    yl = ylim;
 
    % 95th percentile shading
    x_shade = [threshold_95, max(null_dist), max(null_dist), threshold_95];
    y_shade = [0, 0, yl(2), yl(2)];
    fill(x_shade, y_shade, [0.9 0.7 0.7], ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
        'DisplayName', '95th percentile');
 
    % Null mean (the null is a null of per-fly means, so compare means for consistency)
    plot([null_mean null_mean], yl, 'k--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Null mean (%.3f)', null_mean));
 
    plot([obs_mean obs_mean], yl, '-', ...
        'Color', colors.observed, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Observed mean (%.3f)', obs_mean));
 
    plot([0 0], yl, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
        'HandleVisibility','off');
 
    xlabel('Mean OLS slope');
    ylabel('Count');
    title(sprintf(['Population permutation null (slope)\n' ...
        'Null mean = %.3f\np = %.4f (n = %d flies)'], ...
        null_mean, p_val, nFlies));
    legend('Location', 'northwest', 'FontSize', 8);
    box off;
    set(gca, 'TickDir', 'out');
end


% ---
function plot_population_signal_variability(S)

    hold on;

    nFlies = numel(S.fly_IDs_all);

    % Distinct colour per fly - CVD-safe palette
    okabe_ito = [ 
        0.90 0.60 0.00;   % orange
        0.35 0.70 0.90;   % sky blue
        0.00 0.60 0.50;   % bluish green
        0.95 0.90 0.25;   % yellow
        0.00 0.45 0.70;   % blue
        0.80 0.40 0.00;   % vermillion
        0.80 0.60 0.70;   % purple
        0.00 0.00 0.00];  % black
    
    cmap = okabe_ito(mod((1:max(nFlies,7))-1, 8) + 1, :);
    %cmap = lines(max(nFlies,7));

    legend_handles = gobjects(0);

    for f = 1:nFlies

        if f > numel(S.proj_actual_dist) || ...
           isempty(S.proj_actual_dist{f}) || ...
           isempty(S.r_orth_actual_dist{f})
            continue
        end

        x = S.proj_actual_dist{f};
        y = S.r_orth_actual_dist{f};

        valid = isfinite(x) & isfinite(y);

        x = x(valid);
        y = y(valid);

        if numel(x) < 3
            continue
        end

        colour = cmap(f,:);

        % Scatter points
        h = scatter(x, y, 25, ...
            'filled', ...
            'MarkerFaceColor', colour, ...
            'MarkerFaceAlpha', 0.3);

        legend_handles(end+1) = h;

        % OLS fit
        mdl = fitlm(x, y);

        xfit = linspace(min(x), max(x), 100)';
        yfit = mdl.Coefficients.Estimate(1) + ...
               mdl.Coefficients.Estimate(2)*xfit;

        plot(xfit, yfit, ...
            'Color', colour, ...
            'LineWidth', 2);
    end

    xlabel('Response amplitude');
    ylabel('Response variability');

    title('Actual MBON data: response variability vs response amplitude');
    set(gca, 'TickDir', 'out');
    box off;
    axis square;

    set(gca, 'Position', [0.08, 0.12, 0.70, 0.78]);
    lgd = legend(legend_handles, S.fly_IDs_all, 'Interpreter', 'none', ...
        'FontSize', 7);
    lgd.Position = [0.80, 0.15, 0.19, 0.7];
    lgd.Box = 'on';
    

end


% ---
function plot_sn_matched_control(slope_dist_pooled, rho_dist_pooled, ...
                                  slope_obs, rho_obs, ...
                                  p_slope, p_rho, flyID)
% Per-fly matched-sample shot noise control plot
% Two panels: OLS slope distribution (left) and rho (right)
% Pooled ROI1 + ROI2 distribution shown as histogram; actual value in red
% Mirrors layout of plot_permutation_null

    panel_data = { ...
        slope_dist_pooled, slope_obs, p_slope, 'OLS slope'; ...
        rho_dist_pooled, rho_obs, p_rho, 'Spearman rho'};

    for i = 1:2
        subplot(1, 2, i);
        hold on;

        dist    = panel_data{i, 1};
        obs_val = panel_data{i, 2};
        p_val   = panel_data{i, 3};
        stat_lbl = panel_data{i, 4};

        histogram(dist, 50, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', ...
            'DisplayName', 'Matched SN dist (ROI1+ROI2 pooled)');

        yl = ylim;

        % Observed actual value
        plot([obs_val obs_val], yl, 'r-', 'LineWidth', 2.5, ...
            'DisplayName', sprintf('Actual (%.3f)', obs_val));

        % Median of pooled distribution
        med_val = median(dist, 'omitnan');
        plot([med_val med_val], yl, 'k--', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('SN median (%.3f)', med_val));

        xlabel(stat_lbl);
        ylabel('Count');
        title(sprintf('Matched SN null: %s\np_{matched} = %.4f', stat_lbl, p_val), ...
            'Interpreter', 'tex');
        legend('Location', 'northwest', 'FontSize', 8);
        box off;
        set(gca, 'TickDir', 'out');
    end
    box off;
    sgtitle(sprintf('Fly %s: Sample-size matched shot noise control\n(%d pooled values)', ...
        flyID, numel(slope_dist_pooled)), 'Interpreter', 'none');
end


% ---
function plot_sn_matched_population(S)
% Population-level matched shot noise null test
% Monte Carlo null distribution vs observed mean actual slope


    null_dist    = S.pop.sn_matched_slope.null_dist;
    obs_mean     = S.pop.sn_matched_slope.obs_mean;
    null_mean    = S.pop.sn_matched_slope.null_mean;
    p_val        = S.pop.sn_matched_slope.p_pop;
    threshold_95 = S.pop.sn_matched_slope.threshold_95;
    nFlies       = S.pop.sn_matched_slope.nFlies;

    colors.null      = [0.7 0.7 0.7];
    colors.observed  = [0.8 0.2 0.2];

    hold on;

    histogram(null_dist, 60, ...
        'FaceColor', colors.null, 'EdgeColor', 'none', ...
        'DisplayName', 'Matched SN null (ROI1+ROI2 pooled)');

    yl = ylim;

    % 95th percentile shading
    x_shade = [threshold_95, max(null_dist), max(null_dist), threshold_95];
    y_shade = [0, 0, yl(2), yl(2)];
    fill(x_shade, y_shade, [0.9 0.7 0.7], ...
        'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
        'DisplayName', '95th percentile');

    % Null mean (the null is a null of per-fly means, so compare means for consistency)
    plot([null_mean null_mean], yl, 'k--', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Null mean (%.3f)', null_mean));

    % Observed mean actual slope
    plot([obs_mean obs_mean], yl, '-', ...
        'Color', colors.observed, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Observed mean (%.3f)', obs_mean));

    % Zero reference
    plot([0 0], yl, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
        'HandleVisibility', 'off');

    xlabel('Mean OLS slope');
    ylabel('Count');
    title(sprintf(['Population matched-SN null (slope)\n' ...
        'Null mean = %.3f\np = %.4f (n = %d flies)'], ...
        null_mean, p_val, nFlies));
    legend('Location', 'northwest', 'FontSize', 8);
    box off;
    set(gca, 'TickDir', 'out');
end


% ---
function plot_combined_null_violins_perm(S)
% Permutation null violin figure

    [x_lims, actual_slopes, fly_labels, x_pop, x_flies, ...
     col_pop_null, col_fly_null, col_obs] = violin_shared_setup(S);

    pop_dist  = S.pop.slope_perm.null_dist;
    obs_mean  = S.pop.slope_perm.obs_mean;
    fly_dists = S.perm_null_slope_dist;
    p_val     = S.pop.slope_perm.p_pop_slope;

    violin_figure_core(pop_dist, obs_mean, fly_dists, actual_slopes, ...
        x_pop, x_flies, fly_labels, x_lims, ...
        col_pop_null, col_fly_null, col_obs, col_obs, ...
        'Population null dist', 'Per-fly null dist', ...
        sprintf('Permutation null  |  Population p = %.4f', p_val));
end

% ---
function plot_combined_null_violins_sn(S)
% Matched SN control violin figure

    [x_lims, actual_slopes, fly_labels, x_pop, x_flies, ...
     col_pop_null, col_fly_null, col_obs] = violin_shared_setup(S);

    pop_dist  = S.pop.sn_matched_slope.null_dist;
    obs_mean  = S.pop.sn_matched_slope.obs_mean;
    fly_dists = S.slope_dist_sn_matched_pooled_all;
    p_val     = S.pop.sn_matched_slope.p_pop;

    violin_figure_core(pop_dist, obs_mean, fly_dists, actual_slopes, ...
        x_pop, x_flies, fly_labels, x_lims, ...
        col_pop_null, col_fly_null, col_obs, col_obs, ...
        'Population SN null dist', 'Per-fly SN null dist', ...
        sprintf('Matched sample-size SN control  |  Population p = %.4f', p_val));
end

% ---
function violin_figure_core(pop_dist, obs_mean, fly_dists, actual_slopes, ...
        x_pop, x_flies, fly_labels, y_lims, ...
        col_pop_null, col_fly_null, col_pop_obs, col_fly_obs, ...
        leg_pop_null, leg_fly_null, fig_title)
% Shared rendering logic for both violin figures.
% Slope on y-axis; population violin at x=1, per-fly violins at x=3,4,...
% y_lims is the shared slope range computed across both figures.

    nFlies = numel(x_flies);
    hold on;

    % --- Population violin ---
    draw_violin(pop_dist, x_pop, col_pop_null, 0.5);

    % Population observed mean: horizontal line spanning violin
    plot([x_pop - 0.42, x_pop + 0.42], [obs_mean obs_mean], ...
        '-', 'Color', col_pop_obs, 'LineWidth', 2.5);

    % --- Per-fly violins ---
    for f = 1:nFlies
        dist_f = fly_dists{f};
        if isempty(dist_f) || sum(isfinite(dist_f)) < 10; continue; end
        draw_violin(dist_f, x_flies(f), col_fly_null, 0.4);
    end

    % --- Per-fly actual slopes (magenta dots) ---
    for f = 1:nFlies
        if ~isfinite(actual_slopes(f)); continue; end
        scatter(x_flies(f), actual_slopes(f), 70, ...
            col_fly_obs, 'filled', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.95);
    end

    % --- Dividing line between population and per-fly ---
    xline(x_flies(1) - 1, '--', 'Color', [0.82 0.82 0.82], 'LineWidth', 1);

    % --- Zero reference ---
    yline(0, ':', 'Color', [0.60 0.60 0.60], 'LineWidth', 1);

    % --- Axes ---
    x_tick_pos    = [x_pop, x_flies];
    x_tick_labels = [{'Population'}, fly_labels];

    % Major ticks: round y_lims to nearest 0.5, step by 0.05
    y_major_start = floor(y_lims(1) / 0.5) * 0.5;
    y_major_end   = ceil(y_lims(2) / 0.5) * 0.5;
    y_ticks_major = y_major_start:0.05:y_major_end;


    set(gca, ...
        'YLim',       y_lims, ...
        'XLim',       [x_pop - 0.8, x_flies(end) + 0.8], ...
        'XTick',      x_tick_pos, ...
        'YTick',      y_ticks_major, ...
        'XTickLabel', x_tick_labels, ...
        'TickDir',    'out', ...
        'FontSize',   9, ...
        'YMinorTick', 'on');

   % Minor tick spacing: MATLAB sets minor ticks as N subdivisions between majors
    % 0.05 major step / 0.01 minor spacing = 5 minor intervals per major step
    ax = gca;
    ax.YAxis.MinorTickValues = y_major_start:0.01:y_major_end;

    xtickangle(30);
    ylabel('OLS slope');
    title(fig_title, 'FontSize', 10);
    box off;

    % --- Legend (condensed to 2 entries) ---
    col_null_legend = (col_pop_null + col_fly_null) / 2;   % midpoint shade, represents the null family
    h_null = patch(NaN, NaN, col_null_legend, 'EdgeColor', 'none', 'FaceAlpha', 0.65);
    
    h_obs = plot(NaN, NaN, '-o', 'Color', col_pop_obs, ...
        'MarkerFaceColor', col_pop_obs, 'MarkerEdgeColor', 'none', 'LineWidth', 2.5);
    
    legend([h_null, h_obs], {'Null distribution', 'Observed'}, ...
        'Location', 'northeast', 'FontSize', 8);
end


% ---
function [y_lims, actual_slopes, fly_labels, x_pop, x_flies, ...
          col_pop_null, col_fly_null, col_obs] = violin_shared_setup(S)
% Computes shared y-limits, layout, and colours used by both violin figures.

% NB: Apply y-limit - percentile trim to null distribution only. Avoid
% trimming observed values

    nFlies        = numel(S.fly_IDs_all);
    actual_slopes = S.slope_actual_all;
    fly_labels    = arrayfun(@(i) sprintf('Fly%d', i), 1:nFlies, 'UniformOutput', false);

    x_pop   = 1;
    x_flies = (1:nFlies) * 1.5 + 2;

    % Pool null distributions (population + per-fly, both controls)
    null_vals = [ ...
        S.pop.slope_perm.null_dist(:); ...
        S.pop.sn_matched_slope.null_dist(:)];

    for f = 1:nFlies
        d1 = S.perm_null_slope_dist{f};
        d2 = S.slope_dist_sn_matched_pooled_all{f};
        if ~isempty(d1); null_vals = [null_vals; d1(:)]; end
        if ~isempty(d2); null_vals = [null_vals; d2(:)]; end
    end

    null_vals = null_vals(isfinite(null_vals));
    obs_vals  = actual_slopes(isfinite(actual_slopes));

    % Trim null tails, but ensure that observed values are  in range
    null_bounds = prctile(null_vals, [0.1,99.9]);

    y_low  = min([null_bounds(1); obs_vals(:); 0]);
    y_high = max([null_bounds(2); obs_vals(:)]);

    pad    = (y_high - y_low) * 0.05;
    y_lims = [y_low - pad, y_high + pad];

    % Colours
    c = okabe_ito_colors();
    col_pop_null = darken_color(c.blue, 0.55);
    col_fly_null = lighten_color(c.blue, 0.75);
    col_obs      = c.orange;
end


% ---
function draw_violin(data, x_pos, face_col, violin_half_width)
% Draws a single vertical violin at x_pos.
%   data             : vector of values
%   x_pos            : scalar x position
%   face_col         : [r g b] fill colour
%   violin_half_width: half-width of violin at widest point (x-axis units)

    data = data(isfinite(data));
    if numel(data) < 5; return; end

    [f, xi] = ksdensity(data);

    % Normalise so max density reaches violin_half_width
    f_norm = f / max(f) * violin_half_width;

    % Closed polygon: right side (increasing x) then left side (mirrored)
    x_poly = [x_pos + f_norm,  fliplr(x_pos - f_norm)];
    y_poly = [xi,               fliplr(xi)];

    fill(x_poly, y_poly, face_col, ...
        'EdgeColor', face_col * 0.6, ...
        'FaceAlpha', 0.65, ...
        'LineWidth', 0.8);

    % Mean tick
    ctr_val   = mean(data);
    [~, cidx] = min(abs(xi - ctr_val));
    half_w    = f_norm(cidx);
    plot([x_pos - half_w, x_pos + half_w], [ctr_val, ctr_val], ...
        '-', 'Color', face_col * 0.5, 'LineWidth', 2);
end


% --- Helpers ---
function ax = make_narrow_subplot(m, n, p, scale)
    ax = subplot(m, n, p);
    pos = get(ax, 'Position');
    new_width = pos(3) * scale;
    shift = (pos(3) - new_width) / 2;
    set(ax, 'Position', [pos(1)+shift, pos(2), new_width, pos(4)]);
end

function beautify_axis()
    box off; set(gca, 'TickDir', 'out', 'LineWidth', 1);
    axis tight; 
    yl = ylim;
    ylim([yl(1)-abs(yl(1)*0.05), yl(2)+abs(yl(2)*0.05)]);
end

function s = pstars(p)
    if isnan(p), s = 'n.s.';
    elseif p < 0.001, s = '***';
    elseif p < 0.01, s = '**';
    elseif p < 0.05, s = '*';
    else, s = 'n.s.'; 
    end
end

function out = darken_color(rgb, amount)
    % amount in [0,1]; higher = darker
    out = rgb * (1 - amount);
end

function out = lighten_color(rgb, amount)
    % amount in [0,1]; higher = lighter
    out = rgb + (1 - rgb) * amount;
end

%% ----------------------------------------------------------------------------
%                               HELPER FUNCTIONS - ANALYSIS
% -----------------------------------------------------------------------------


% --- TLS REGRESSION ---
function [v_tls, proj_sign, r_orth] = tls_projection(x, y, label, silent, input_lambda)
   
    if nargin < 3 
        label = ''; 
    end
    if nargin < 4
        silent = false;
    end
   
  
    X_raw = [x(:), y(:)]; % ensure data points stay paired before 'cleaning'
    valid = all(isfinite(X_raw), 2); % logical column vector marking rows where bith x and y are valid
    x_clean = X_raw(valid, 1); % extract only valid rows 
    y_clean = X_raw(valid, 2);
    mu = mean([x_clean, y_clean], 1);

    % Ensure Lambda is provided
    if nargin < 5 || isempty(input_lambda)
        error('Execution Halted: Lambda not specified for "%s". Within-stimulus noise ratio is required.', label);
    end
    
    lambda = input_lambda;

    % --- DEMING REGRESSION CALCULATION ---
    % Calculate covariance
    s_xy = cov(x_clean, y_clean); 
    s_xx = s_xy(1,1);
    s_yy = s_xy(2,2);
    s_xy_val = s_xy(1,2);
    
    % Deming Slope Formula
    num = s_yy - lambda * s_xx + sqrt((s_yy - lambda * s_xx)^2 + 4 * lambda * s_xy_val^2);
    den = 2 * s_xy_val;
    m = num / den;     % Slope
    
    % Handle perfect vertical/horizontal cases
    if abs(s_xy_val) < 1e-9
        if s_xx > s_yy
            m = 0; 
        else
            m = Inf; 
        end
    end

    % Vector logic, ensure v_tls is unit vector a^2 + b^2 = 1 
    % (v_tls = [a; b], where a = x-component (MBONap)and b = y-component(MBONav))
    v_tls = [1; m];
    v_tls = v_tls / norm(v_tls);

    % Enforce positive direction (pointing to top-right)
    % Using proj_sign as proxy for signal strength, so large projection
    % value = high response amplitude (not mathematical requirement, just for specific dataset)
    if v_tls(1) < 0
    v_tls = -v_tls;
    end

    v_tls = v_tls(:); % Enforce column vector, regardles on how constructed
    
    % --- Calculate Residuals ---

    % Center data
        % Each point shifted so that centroid is at origin, 
        % Each point is expressed as deviation from mean response pair
    Xc = [x_clean - mu(1), y_clean - mu(2)]; 
    
    % Projection
        % Scalar projection of each centred data point onto TLS unit vector
        % i.e., how far along TLS axis that point sits, measured from centroid
    proj_from_mean = Xc * v_tls; 

    % Residuals
        % Pythagorean theorem in vector spaces — any vector in 2D can be uniquely decomposed into two 
        % perpendicular components, and those components add up exactly to the original vector 
        % Xc = component along axis (proj_from_mean * v_tls') + component perpendiculat to axis (orth_vectors)
    orth_vectors = Xc - (proj_from_mean * v_tls'); 
    
    % Unsigned
    r_orth = sqrt(sum(orth_vectors.^2, 2)); % Euclidean length of each perpendicular vector = unsigned residual

    % Signed (|r_signed| = r_orth) - not currently used
    %v_perp = [-v_tls(2); v_tls(1)]; % v_perp is the unit vector perpendicular to TLS axis.
    %                                  rotating [a;b] by 90deg gives [-b;a]. ([a;b] . [-b;a] = 0 
    %                                  (confirms perpendicularity))
    %r_signed = Xc * v_perp; % Projecting directly onto v_perp, rather than TLS axis, 
    %                          gives a signed scalar - distance from TLS
    %                          axis, plus which side point falls on.

    % Shift "Drive" so that proj_sign reflects absolute response magnitude, 
    % not deviation from the mean. Important for the interpretation and ordering 
    % of trials along the x-axis (from weakest to strongest response, not from 
    % most-negative-deviation to most-positive-deviation. 
    % NB: 'Spearman rho' is invariant and 'slope' is invariant to a shift
    % in x -> only intercept changes.
   
    min_val = min(proj_from_mean);
    proj_sign = proj_from_mean - min_val; 
    
    if ~silent
    fprintf('Regression (%s): Lambda=%.3f, Slope=%.3f\n', label, lambda, m);
    end

    
    % Compare Deming 'correction' to raw TLS via SVD
    % The angle between Deming and raw SVD axes is key diagnostic:
    % Small angle (<5°) — lambda is close to 1, Deming barely differs from TLS, correction is minor
    % Moderate angle (5-20°) — meaningful correction, Deming is doing useful work
    % Large angle (>20°) — strong correction, verify lambda is reasonable and check both axes make biological sense
    % Very large angle (>35°) — Scrutinising carefully

    [~,~,V] = svd(Xc, 'econ');
    v_svd = V(:,1);
    if v_svd(1) < 0 
        v_svd = -v_svd; 
    end

    if ~silent
    fprintf('Angle vs raw TLS (deg): %.4f\n', acosd(min(1, abs(dot(v_tls, v_svd)))));
    fprintf('Deming axis:        [%.4f, %.4f]\n', v_tls(1), v_tls(2));
    fprintf('Raw TLS axis:       [%.4f, %.4f]\n', v_svd(1), v_svd(2));
    fprintf('Corr(proj, r_orth): %.6f\n',         corr(proj_sign, r_orth));
    fprintf('Corr(proj^2, r_orth^2): %.6f\n',     corr(proj_sign.^2, r_orth.^2));
    end
end


    
% --- PAIRED PERMUTATION TEST ---
function pval = paired_perm(x, y, statFun, nPerm)
    % Keep only finite values
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    
    % Convert to column vectors to avoid shape issues
    x = x(:);
    y = y(:);
    
    % Difference vector
    d = x - y;
    
    % Compute observed statistic
    obsStat = statFun(d);
    if ~isscalar(obsStat)
        error('statFun must return a scalar for the observed statistic.');
    end
    
    n = numel(d);
    
    % Preallocate permutation stats
    permStats = zeros(nPerm,1);
    
    for k = 1:nPerm
        % Random ±1 signs
        signs = randi([0 1], n, 1)*2 - 1;
        
        % Ensure element-wise multiplication produces a column vector
        temp = d .* signs;
        
        % Compute statistic for this permutation
        testVal = statFun(temp);
        
        % Defensive check: must be scalar
        if ~isscalar(testVal)
            error('statFun must return a scalar for permutation %d', k);
        end
        
        permStats(k) = testVal;
    end
    
    % One-sided p-value
    pval = mean(permStats >= obsStat);
end