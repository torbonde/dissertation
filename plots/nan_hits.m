function X = nan_hits(Y, L, U)
X = Y;
I = [];
[~, N] = size(Y);
for i = 1:N-1
    I1 = find(U < X(:,i) | X(:,i) < L);
    I = union(I,I1);
    X(I, i) = NaN;
end
X(I,end) = NaN;