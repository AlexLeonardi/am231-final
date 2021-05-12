function X_grey = double_to_image2( X_double )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[m, n] = size(X_double);
X_grey = zeros(m, n);
for i = 1:m
    for j = 1:n
        X_grey(i, j) = round(X_double(i, j));
    end
end

X_grey = uint8(X_grey);

end

