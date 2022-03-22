function [CI_low, CI_high] = confidence_interval(X, dim, confidence)
    if(nargin < 2 || isempty(dim))
        dim = 1;
    end
    if(nargin < 3 || isempty(confidence))
        confidence = .025;
    end
    X = sort(X, dim); n = size(X, dim); side = round(n * confidence);
    switch dim
        case 1
            CI_low = X(side+1,:,:); CI_high = X(n - side,:,:);
        case 2
            CI_low = X(:,side+1,:); CI_high = X(:,n - side,:);
        case 3
            CI_low = X(:,:,side+1); CI_high = X(:,:, n-side);
        otherwise
            fprintf('bad dimension!')
    end
end