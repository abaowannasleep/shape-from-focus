function ssdErr = ssd(fm, focus, A, zi, s)
[M, N, P] = size(fm);
ssdErr =  zeros(M, N);
    for p = 1:P
        err = fm(:,:,p) - A.*exp(- (focus(p)- zi).^2 ./ (2*s.^2) );
        ssdErr = ssdErr + err.*err;  
    end 
    ssdErr = ssdErr/P;
end
