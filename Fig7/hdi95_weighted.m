function hdi = hdi95_weighted(samples, weights)
    samples = samples(:);
    weights = weights(:) / sum(weights);
    
    [samples, idx] = sort(samples);
    weights = weights(idx);
    
    cdf = cumsum(weights);
    
    credMass = 0.95;
    minWidth = Inf;
    hdi = [NaN NaN];
    
    for i = 1:length(samples)
        j = find(cdf >= cdf(i) + credMass, 1);
        if ~isempty(j)
            width = samples(j) - samples(i);
            if width < minWidth
                minWidth = width;
                hdi = [samples(i), samples(j)];
            end
        end
    end
end
