function [CI_low, CI_high] = confidence_interval(X, confidence) %goes by columns
if(nargin < 2)
    confidence = .025;
end
X = sort(X, 1); s= size(X); n = s(1); side = round(n * confidence);
CI_low = X(side+1,:,:); CI_high = X(n - side,:,:);
end