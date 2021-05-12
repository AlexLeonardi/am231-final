function  vec_return = A_operator_T( Omega, mn, vec)
    vec_return = zeros(mn, 1);
    for i = 1:length(Omega)
        vec_return(Omega(i)) = vec(i);
    end    
end

