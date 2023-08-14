% codes ref https://www.mathworks.com/matlabcentral/fileexchange/55103-shape-from-focus
function ks_r = kstest(fm, focus, I, A, zi, s)
    fmax = focus(I);
    [M, N, P] = size(fm);
    err = zeros(M, N);
    for p = 1:P
        err = err + abs( fm(:,:,p) - A.*exp(- (focus(p)- zi).^2 ./ (2*s.^2) ) );
    end 
    ks_r = -real(20*log10(P*fmax./err));
end