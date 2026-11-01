function weights = drawWeights_PWgivenNtheta(MBmodel)
    W                    = load('data/W_PN_KC.mat').W;
    PW_given_theta_and_n = load('data/PW_given_theta_and_n.mat').PW_given_theta_and_n;
    weights = zeros(MBmodel.nPNs, MBmodel.nKCs);
    for k = 1:MBmodel.nKCs
        conn     = MBmodel.connectivity{k};
        N        = numel(conn);
        ThetaInd = round((MBmodel.theta(k)-0.01)/0.1 + 1);
        ThetaInd(ThetaInd<1) = 1;                       
        pW       = squeeze(PW_given_theta_and_n(N-1, ThetaInd, :));
        w        = randsample(W, N, 'true', pW);        % one draw per claw (pW is constant within a KC)
        weights(:,k) = accumarray(conn(:), w(:), [MBmodel.nPNs,1]);  % accumulate duplicates
    end
end