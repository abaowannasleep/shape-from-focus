function [d1, err_list] = gradient_descent(y_obs, d, mask, alpha)
n = 1;
labmda = 0.2;
labmda2 = 1;
err_list = zeros([1, 200]);
g = 2*( mask.*( (d-y_obs)+labmda*pair_cost(d) ) + (1-mask).*(labmda2*pair_cost(d)) );  
% g = 2 * (mask.*(mask .* d - y_obs) + labmda*pair_cost(d)) ;
d1 = d - alpha .* g;
err = sum( sum((mask.*(d - d1)).^2) );
err_list(n) = err;
while( n < 200)
    d = d1;
%     g = 2 * (mask.*(mask .* d - y_obs) + labmda*pair_cost(d)) ;
    g = 2*( mask.*( (d-y_obs)+labmda*pair_cost(d) ) + (1-mask).*((1-labmda)*pair_cost(d)) ); 
    d1 = d - alpha .* g;
    n = n + 1;
    err = sum( sum((mask.*(d - d1)).^2) );
    err_list(n) = err;
end
end

function pc = pair_cost(d)
kernel = [[-0.5 -1, -0.5];[-1, 6, -1];[-0.5, -1, -0.5 ]];
pc = imfilter(d, kernel, 'replicate');
end