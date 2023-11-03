function I=allinfocus(d,xx)
 [w, h] = size(d);
 I = zeros(w,h);
    for i=1:w
        for j=1:h
            I(i,j)=xx(i,j,d(i,j));
        end
    end
end