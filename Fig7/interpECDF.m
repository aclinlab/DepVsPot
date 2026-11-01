function [y, bins] = interpECDF(data, bins)
% [y, bins] = interpECDF(data, bins)
% This function takes an empirical cumulative distribution function (ECDF)
% and interpolates the x-axis onto the specified points
% Inputs:
% data: any data that the function ecdf can deal with. The data for which
% you want to draw the CDF.
% bins: any data that can serve as the input Xq in the function interp1.
% This would be the x-values that you want to use to interpolate the ECDF
% Returns the same bins used as input

if range(data(isfinite(data)))==0 % all elements are equal - ignoring NaN and Inf
    % this is required because when all data are equal, the ECDF is a
    % vertical line, and it's impossible to do interp1 on a vertical line
    value = unique(data(isfinite(data)));
    y = zeros(length(bins),1);
    y(bins<value) = 0;
    y(bins>=value) = 1;
else
    [f,x] = ecdf(data);
    y = interp1(x(2:end), f(2:end), bins, 'previous');
end

end