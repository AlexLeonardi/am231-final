function  vec_return = A_operator( Omega, matrix_vec)
    vec_return = zeros(length(Omega), 1);
    for i = 1:length(Omega)
        vec_return(i) = matrix_vec(Omega(i));
    end    
end

