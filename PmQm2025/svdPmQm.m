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
alpha = 0.25;
f     = @(x, y) (x.^2+y.^2).^alpha;

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
        yx = @(x) -x;
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

% L2 SVD, H1-L2 and L2-H1 SVD
[U00, S00, V00] = svd(C);
[U10, S10, V10] = svd(C10);
[U01, S01, V01] = svd(C01);
s00             = diag(S00);
s10             = diag(S10);
s01             = diag(S01);

% Truncate
errL2        = zeros(size(C, 1)+1, 1);
errL2L2      = zeros(size(C, 1)+1, 1);
errH1L2      = zeros(size(C, 1)+1, 1);
errL2H1      = zeros(size(C, 1)+1, 1);
errH1L2_L2H1 = zeros(size(C, 1)+1, 1);
errsvs       = zeros(size(C, 1)+1, 1);

n  = 2*K+1;
for r  = 0:1:n
    % L2 SVD error
    errL2(r+1)      = norm(s00(r+1:end));

    % H1L2 and L2H1 errors
    errsvs(r+1) = norm([s10(r+1:end); s01(r+1:end)]);

    % Scaling
    k    = -K:K;
    S    = (1+k.^2).^(-1/2);
    S    = diag(S);
    Sinv = (1+k.^2).^(1/2);
    Sinv = diag(Sinv);

    % L2-L2 projection
    U   = U00(:, 1:r);
    V   = V00(:, 1:r);
    P   = kron(V*V', U*U');
    c   = C(:);
    D00 = P*c;
    D00 = reshape(D00, n, n);

    % H1-L2 projection
    U   = U10(:, 1:r);
    V   = V10(:, 1:r);
    P10 = kron(V*V', U*U');
    c   = Sinv*C;
    c   = c(:);
    D10 = P10*c;
    D10 = reshape(D10, n, n);
    D10 = S*D10;

    % L2-H1 projection
    U   = U01(:, 1:r);
    V   = V01(:, 1:r);
    P01 = kron(V*V', U*U');
    c   = C*Sinv;
    c   = c(:);
    D01 = P01*c;
    D01 = reshape(D01, n, n);
    D01 = D01*S;

    % H1-L2*L2-H1 projection
    c      = Sinv*D01;
    c      = c(:);
    D10_01 = P10*c;
    D10_01 = reshape(D10_01, n, n);
    D10_01 = S*D10_01;

    % Compute H1-error
    delta   = zeros(n, n);
    for k = 1:n
        for m = 1:n
            delta(k, m)   = (C(k, m)-D00(k, m))*...
                             sqrt(1+(-K+k-1)^2+(-K+m-1)^2);

        end
    end
    errL2L2(r+1) = norm(delta, 'fro');

    for k = 1:n
        for m = 1:n
            delta(k, m)   = (C(k, m)-D10(k, m))*...
                             sqrt(1+(-K+k-1)^2+(-K+m-1)^2);

        end
    end
    errH1L2(r+1) = norm(delta, 'fro');

    for k = 1:n
        for m = 1:n
            delta(k, m)   = (C(k, m)-D01(k, m))*...
                             sqrt(1+(-K+k-1)^2+(-K+m-1)^2);

        end
    end
    errL2H1(r+1) = norm(delta, 'fro');

    for k = 1:n
        for m = 1:n
            delta(k, m)   = (C(k, m)-D10_01(k, m))*...
                             sqrt(1+(-K+k-1)^2+(-K+m-1)^2);

        end
    end
    errH1L2_L2H1(r+1) = norm(delta, 'fro');
end

% Write error vectors
dlmwrite('errL2.dat', errL2)
dlmwrite('errL2L2.dat', errL2L2)
dlmwrite('errH1L2.dat', errH1L2)
dlmwrite('errL2H1.dat', errL2H1)
dlmwrite('errH1L2_L2H1.dat', errH1L2_L2H1)
dlmwrite('errsvs.dat', errsvs)

fprintf('Output written, programm done\n')
