function [I,d] = computeST2D(img)
    [X,Y,Z]=size(img);
    
    t11=zeros(X,Y,Z);
    t12=zeros(X,Y,Z);
    t22=zeros(X,Y,Z);
    for z=1:Z
        [t11(:,:,z), t12(:,:,z), t22(:,:,z)] = compute_structure_tensor2d(img(:,:,z));
    end
    
    % T=t22;
    % T=t11+t12+t22;
    T=t11+t22;
    [~,d]=max(T,[],3);
    I=allinfocus(d,img);
end