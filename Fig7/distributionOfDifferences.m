function result = distributionOfDifferences(x,xw,y,yw,k)
% given distributions x and y, weighted by xw and yx, take k random pairs
% and get the difference between x and y for each pair

index1 = randsample(size(x,1),k,'true',xw);
index2 = randsample(size(y,1),k,'true',yw);
result = x(index1)-y(index2);

end