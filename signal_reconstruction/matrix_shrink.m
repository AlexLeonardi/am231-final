function return_matrix = matrix_shrink( U, S, V, nu)
    s = diag(S);
    s_hat = vector_shrink(s, nu);    
    S_hat = diag(s_hat);
    
    return_matrix = U*S_hat*V';
        
end

