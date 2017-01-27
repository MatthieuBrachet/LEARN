function P = phillips(kx, ky, windDir, windSpeed, A, g)
    k_sq = kx^2 + ky^2;
    L = windSpeed^2 / g;
    k = [kx, ky] / sqrt(k_sq);
    wk = k(1) * windDir(1) + k(2) * windDir(2);
    P = A / k_sq^2 * exp(-1.0 / (k_sq * L^2)) * wk^2;
end