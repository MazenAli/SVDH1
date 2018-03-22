% Zero precision for linear dependence
eps = 1e-09;

% Domain
xmin  = -pi;
xmax  = pi;
ymin  = -pi;
ymax  = pi;

% Fourier basis
psik  = @(x, k) exp(-sqrt(-1)*k*x);

% Truncation parameters
K  = 30;

% Test function
alpha = 0.6;
c     = 0.5;
f     = @(x, y) abs(c*x+y).^alpha;

% Compute coefficients
C   = sqrt(-1)*ones(2*K+1, 2*K+1);
C10 = C;
C01 = C;

fprintf('\n\n')
for k = -K:K
    for m=-K:K
        fprintf('Calculating k=%d, m=%d\n', k, m)

        scale = 2*pi;
        fun   = @(x, y)(f(x,y).*psik(x, k).*psik(y, m)/scale);

        % Split around singularity
        yx = @(x) -c*x;
        val = integral2(fun, xmin, xmax, ymin, yx);
        val = val+integral2(fun, xmin, xmax, yx, ymax);

        C(k+K+1, m+K+1)   = val;
        C10(k+K+1, m+K+1) = val*sqrt(1+k^2);
        C01(k+K+1, m+K+1) = val*sqrt(1+m^2);
    end
end

% Test truncation quality
p = 0.2;
M = ceil(K*p);
D = abs(C);
D(M+1:2*K+1-M, M+1:2*K+1-M) = zeros(2*(K-M)+1, 2*(K-M)+1);

rel = norm(D, 'fro')/norm(C, 'fro');
fprintf('The relative error is %f\n', rel)

% SVD

% L2 SVD
s = svd(C);

% H1-L2 and L2-H1 SVD
[U10, S10, V10] = svd(C10);
[U01, S01, V01] = svd(C01);
s10             = diag(S10);
s01             = diag(S01);

% Sort singular values
scup            = [s10 ones(length(s10), 1); s01 2*ones(length(s01), 1)];
scup            = sortrows(scup, -1);

err1            = zeros(size(C, 1)+1, 1);
err2            = zeros(length(scup)+1, 2);
err3            = zeros(size(C, 1)+1, 1);

% Iterate over r=r1+r2;
n  = 2*K+1;
for r  = 0:1:n
    % L2 SVD error
    err1(r+1)            = norm(s(r+1:end));

    % H1-L2 SVD
    S10_                   = S10;
    S10_(r+1:end, r+1:end) = zeros(size(S10, 1)-r, size(S10, 1)-r);
    D10                    = U10*S10_*V10';

    % Determine r1, r2
    id                   = scup(1:r, 2);
    r1                   = sum(id == 1);
    r2                   = sum(id == 2);
    if r>0
        r1 = max(r1, 1);
        r2 = max(r2, 1);
    end

    % Compute scalings
    Dk = -K:K;
    S1 = Dk.^2;
    S2 = (1+Dk.^2).^(-1/2);

    S1 = diag(S1);
    S2 = diag(S2);

    % Compute V, W
    X1 = U10(:, 1:r1);
    Y1 = V01(:, 1:r2);

    % Scale
    X1 = S2*X1;
    Y1 = S2*Y1;

    % Left and right spaces
    X  = X1;
    Y  = Y1;

    % Compute RHS
    RHS = kron(Y', X'+X'*S1)+kron(Y'*S1, X');
    b   = RHS*C(:);

    % Compute system matrix
    A  = RHS*kron(Y, X);
    fprintf('Condition basis r1=%d, r2=%d: kappa = %f\n', r1, r2, cond(A))

    % Solve and transform
    c = [];
    if r>0
        c  = A\b;
    end
    c  = reshape(c, r1, r2);

    % Compute Fourier coefficients
    D  = X*c*Y';

    % Compute H1-error
    delta   = zeros(n, n);
    delta10 = zeros(n, n);
    for k = 1:n
        for m = 1:n
            delta(k, m)   = (C(k, m)-D(k, m))*...
                             sqrt(1+(-K+k-1)^2+(-K+m-1)^2);
            delta10(k, m) = (C(k, m)-D10(k, m)/sqrt(1+(-K+k-1)^2))*...
                             sqrt(1+(-K+k-1)^2+(-K+m-1)^2);

        end
    end

    err2(r+1, 1) = min(r1, r2);
    err2(r+1, 2) = norm(delta,   'fro');
    err3(r+1)    = norm(delta10, 'fro');
end

% Write error vectors
dlmwrite('errL2.dat', err1)
dlmwrite('errH1opt.dat', err2)
dlmwrite('errH1L2.dat', err3)

fprintf('Output written, programm done\n')
