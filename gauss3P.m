function [I, u, s, A] = gauss3P(x, Y)
% Closed-form solution for Gaussian interpolation using 3 points
% Guss func： y = A*e^[- (x-u)^2 / (2*s^2)]
% Internal parameter:
STEP = 2;
[M,N,P] = size(Y);
[~, I] = max(Y,[ ], 3);
[IN,IM] = meshgrid(1:N,1:M);
Ic = I(:);
Ic(Ic<=STEP)=STEP+1;
Ic(Ic>=P-STEP)=P-STEP;
% （maxPos-STEP，maxPos+STEP）
Index1 = sub2ind([M,N,P], IM(:), IN(:), Ic-STEP);
Index2 = sub2ind([M,N,P], IM(:), IN(:), Ic);
Index3 = sub2ind([M,N,P], IM(:), IN(:), Ic+STEP);
x1 = reshape(x(Ic(:)-STEP),M,N);
x2 = reshape(x(Ic(:)),M,N);
x3 = reshape(x(Ic(:)+STEP),M,N);
y1 = reshape(log(Y(Index1)),M,N);
y2 = reshape(log(Y(Index2)),M,N);
y3 = reshape(log(Y(Index3)),M,N);

c = ( (y1-y2).*(x2-x3)-(y2-y3).*(x1-x2) )./...
    ( (x1.^2-x2.^2).*(x2-x3)-(x2.^2-x3.^2).*(x1-x2) );
b = ( (y2-y3)-c.*(x2-x3).*(x2+x3) )./(x2-x3);
a = y1 - b.*x1 - c.*x1.^2;
s = sqrt(-1./(2*c));
u = b.*s.^2;
A = exp(a + u.^2./(2*s.^2));
u(u>max(x)) = max(x);
u(u<min(x)) = min(x);
end
