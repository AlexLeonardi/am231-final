close all
clear

fig_pos_old = [10, 455, 350, 350; 360, 455, 350, 350; 710, 455, 350, 350; 1060, 455, 350, 350; 10, 30, 350, 350; 360, 30, 350, 350; 710, 30, 350, 350; 1060, 30, 350, 350];
fig_pos = [1 498 173 130; 175 498 173 130; 349 498 173 130; 523 498 173 130; 697 498 173 130; 871 498 173 130; 1045 498 173 130];

figure1 = figure('Name', 'Imported Image', 'Renderer', 'painters', 'Position', fig_pos(1,:));
figure2 = figure('Name', 'Grey scale check', 'Renderer', 'painters', 'Position', fig_pos(2,:));
figure3 = figure('Name', 'Truncated Rank', 'Renderer', 'painters', 'Position', fig_pos(3,:));
figure4 = figure('Name', 'Noised Truncated Rank', 'Renderer', 'painters', 'Position', fig_pos(4,:));
figure5 = figure('Name', 'convergence check', 'Renderer', 'painters', 'Position', fig_pos(5,:));
figure6 = figure('Name', 'convergence check2', 'Renderer', 'painters', 'Position', fig_pos(6,:));

% X_rgb = imread('uob_logo_large.png');
X_rgb = imread('uob_logo_small.png');
X_grey_uint8 = rgb2gray(X_rgb);
figure(figure1.Number)
imshow(X_rgb)
X_grey = double(X_grey_uint8);
s_image = svd(X_grey);

X_grey_check = double_to_image2(X_grey);
figure(figure2.Number)
imshow(X_grey_check)

X00 = X_grey;
[U_00,S_00,V_00]=svd(X_grey, 'econ');
[m, n] = size(X00);
r00 = rank(X00);
r = ceil(r00/6);
S_0 = S_00; %% truncate the rank
S_0(r+1:n, r+1:n) = 0;
X0 = U_00* S_0* V_00';
X0_check = double_to_image2(X0);
figure(figure3.Number)
imshow(X0_check)

X0_vec = X0(:);
rand_mat = randi([0, 1], [m, n]);
Omega = find(rand_mat);
M_raw = zeros(length(Omega), 1);
X = zeros(m*n, 1);
count = 1;
for i = Omega'
    M_raw(count) = X0_vec(i);
    X(i) = X0_vec(i);
    count = count + 1;
end
X = reshape(X, [m, n]);
[U,S,V] = svd(X, 'econ');
X_check1 = double_to_image2(X);
figure(figure4.Number)
imshow(X_check1)
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
    
    Xk_check1 = double_to_image(abs(X_k), S(1, 1));
    figure(figure5.Number)
    imshow(Xk_check1)

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
%                 [U_k1, S_k1, V_k1] = svd(X_k, 'econ');
%                 cond2 = norm(U_k1*V_k1' + reshape(G_func_Omega(Omega, M, X_k, b_k)/mu, [m, n]) );
%                 if cond2 < gtol
%                    cond = 1;
%                 end
            cond = 1;
        end

    end
    [breg_it, j, i, rank(X_k)]
    mu = mu*eta;
    svd_last = svd(X_k);
    Xk_check2 = double_to_image(abs(X_k), S(1, 1));
    figure(figure6.Number)
    imshow(Xk_check2)
    norm(X-abs(X_k), 'fro')
    pause(2)
end
norm(X-abs(X_k), 'fro')
Xk_check3 = double_to_image(abs(X_k), S(1, 1));
figure(figure6.Number)
imshow(Xk_check3)

