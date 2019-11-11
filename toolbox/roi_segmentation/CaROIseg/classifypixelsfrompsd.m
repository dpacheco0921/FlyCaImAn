function P = classifypixelsfrompsd(X, pixels2sel)
if ndims(X) > 2
    sizX = size(X);
    X = reshape(X, prod(sizX(1:end-1)), sizX(end));
end
P.active_pixels = zeros(prod(sizX(1:end-1)), 1);
X = bsxfun(@minus, X, mean(X, 2));     % center
X = bsxfun(@times, X, 1./sqrt(mean(X.^2, 2)));
[L,Cx] = kmeans_pp(X(pixels2sel, :)',2);
[~,ind] = min(sum(Cx(max(1,end-49):end,:),1));
if ~isempty(pixels2sel)
   P.active_pixels(pixels2sel) = 1;
   P.active_pixels(pixels2sel(L==ind)) = 2;
else
   P.active_pixels = (L==ind);
end
P.centroids = Cx;
end