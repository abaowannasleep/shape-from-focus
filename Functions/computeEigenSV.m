function [I,d] = computeEigenSV(img)
    [X,Y,Z]=size(img);
    FV=zeros(X,Y,Z);
    d=zeros(X,Y);
    
    WS=9;   % window size
    r=floor(WS/2);  % radius
    N=6;    % number of heighest eigenvalues considered
    
    for z=1:Z
        for x=(r+1):X-r
            for y=(r+1):Y-r
                W=img(x-r:x+r,y-r:y+r,z);   % window
                s=svd(W);
                FV(x,y,z)=sum(s(1:N));
            end
        end
    end
    [~,d]=max(FV,[],3);
    I=allinfocus(d,img);
end
