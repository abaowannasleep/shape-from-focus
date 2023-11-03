function [I,d] = computeST3D(img)

[gx,gy,gz]=gradient(img);

gxx=gx.*gx;
gxy=gx.*gy;
gxz=gx.*gz;
gyy=gy.*gy;
gyz=gy.*gz;
gzz=gz.*gz;

[mu3,mu2,mu1,~,~,~]=eig3volume(gxx, gxy, gxz, gyy, gyz, gzz);

T=abs(mu3)+abs(mu2)+abs(mu1);
% T=abs(mu3)+abs(mu2);
% T=abs(mu3)+abs(mu1);
%  T=abs(mu2)+abs(mu1);

[~,d]=max(T,[],3);
I=allinfocus(d,img);


    
%     sigma=0.6;
    
    %     I = imgaussian(I,sigma);
    
    % [Gx,Gy,Gz]=imgradientxyz(I,'central');
    % [Gxx,Gxy,Gxz]=imgradientxyz(Gx,'central');
    % [Gyx,Gyy,Gyz]=imgradientxyz(Gy,'central');
    % [Gzx,Gzy,Gzz]=imgradientxyz(Gz,'central');
%     
%     [ux,uy,uz] = ModifiedLaplacian(I);
%     
%     Fx2 = imgaussian(ux,sigma);
%     Fy2 = imgaussian(uy,sigma);
%     Fz2 = imgaussian(uz,sigma);
%     FxFy = imgaussian(ux.*uy,sigma);
%     FxFz = imgaussian(ux.*uz,sigma);
%     FyFz = imgaussian(uy.*uz,sigma);
    
    % [Gxx,Gxy,Gxz]=imgradientxyz(Gx);
    % [Gyx,Gyy,Gyz]=imgradientxyz(Gy);
    % [Gzx,Gzy,Gzz]=imgradientxyz(Gz);
    %
    % Gxx = imgaussian(Gx.*Gx,sigma);
    % Gxy = imgaussian(Gx.*Gy,sigma);
    % Gxz = imgaussian(Gx.*Gz,sigma);
    % Gyy = imgaussian(Gy.*Gy,sigma);
    % Gyz = imgaussian(Gy.*Gz,sigma);
    % Gzz = imgaussian(Gz.*Gz,sigma);
    
%     [mu31,mu21,mu11,~,~,~]=eig3volume(Fx2, Fy2, Fz2, FxFy, FxFz, FyFz);
%     
%     stick_tf = abs(mu11) - abs(mu21);
%     plate_tf = abs(mu21) - abs(mu31);
%     ball_tf = abs(mu31);
end