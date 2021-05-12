close all
clear

weeds_num = 10000  
for weeds_it = 239:(weeds_num-1)
    weeds_it
    weed_in = readmatrix(['../deep_weeds_csv_raw/weed' num2str(weeds_it) '.csv']);

    X00 = weed_in(2:size(weed_in, 1),2:size(weed_in, 2) );
    [U_00,S_00,V_00]=svd(X00, 'econ');
    [m, n] = size(X00);
    r00 = rank(X00);
    r = min([r00, ceil(sqrt(m*n  )/6)]);
    S_0 = S_00; %% truncate the rank
    S_0(r+1:n, r+1:n) = 0;
    X0 = U_00* S_0* V_00';

    X0_vec = X0(:);
    rand_mat = randi([0, 1], [m, n]);
    Omega = find(rand_mat);
    M_raw = zeros(length(Omega), 1);
    x = zeros(length(Omega), 1);
    y = zeros(length(Omega), 1);
    v = zeros(length(Omega), 1);
    X = zeros(m*n, 1);
    count = 1;
    for i = Omega'
        M_raw(count) = X0_vec(i);
        X(i) = X0_vec(i);
        count = count + 1;
        x(count) = fix((i-1)/m) + 1;
        y(count) = mod((i-1),m) + 1;
        v(count) = X0_vec(i);
    end
    X = reshape(X, [m, n]);
    [U,S,V] = svd(X, 'econ');

    X_out = uint8(round(X));
    writematrix(X_out, ['../deep_weeds_csv_50pcentMasked/weed' num2str(weeds_it) '.csv']);
    
    [xq,yq] = meshgrid(1:1:n, 1:1:m);
    vq = griddata(x,y,v,xq,yq);
    X_out_grid = uint8(round(max(vq,0)));
    writematrix(X_out_grid, ['../deep_weeds_csv_50pcentMasked_Int/weed' num2str(weeds_it) '_Int.csv'])

    X = X/S(1, 1);
    M = M_raw/S(1, 1);

    b = A_operator(Omega, X(:));
    b_k = b;

    tau = 1;
    mu_low = 1e-8;
    eta = 0.25;

    I_m = 500;

    cond = 1;
    xtol = 1e-10;
    gtol = 1e-4;

    mu = eta*(norm(b_k)^2);

    X_k = zeros(size(X));
    j = 0;
    breg_it = 0;
    while mu > mu_low
        j = j + 1;
        cond = 0;
        nu = mu*tau;
        i = 0;
        while cond == 0 && i < I_m
            i = i + 1;
            X_last = X_k;
            X_025 = X_k(:) - tau *A_operator_T(Omega, m*n, (A_operator(Omega, X_k(:)) - b_k) ) ;
            X_05 = reshape(X_025, [m, n]); %% step 1

            [U_k, S_k, V_k] = svd(X_05, 'econ');

            X_k = matrix_shrink(U_k, S_k, V_k, nu); %% step two

            cond1 = norm(X_k - X_last, 'fro')/max([1, norm(X_last, 'fro') ]);

            if cond1 < xtol
                cond = 1;
            end

        end
%         [breg_it, j, i, rank(X_k)]
        mu = mu*eta;
    end    
    max_check = max(max(X_k*S(1, 1)));
    if max_check > 256
        X_out = uint8(round(X_k*S(1, 1)/max_check * 256));
    else
        X_out = uint8(round(X_k*S(1, 1)));
    end
    
    writematrix(X_out, ['../deep_weeds_csv_50pcentMasked_RR/weed' num2str(weeds_it) '_RR.csv']);

end