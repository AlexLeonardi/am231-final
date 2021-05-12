function vec_hat = vector_shrink( vec , nu)
vec_hat = zeros(length(vec), 1);

for i = 1:length(vec)
    if (vec(i) - nu > 0)
        vec_hat(i) = vec(i) - nu;
    end
    
end


end

