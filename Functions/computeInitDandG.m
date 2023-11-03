function [depth,G1,G2,G3,G4,G5,G6,G7,G8,G9,G10,G11,G12,G13] = computeInitDandG(img)
    
    G1=round(mean(img,3)); % mean intensity in ImgSeq along z direction
    G2=setRange(var(img,0,3));  % variance in ImgSeq along z direction
    
    %% Compute Focus measure Volume and Depth
    %     FV=SML(img);
    FV=GLV(img);
    %     FV=TEN(img);
    [maxF,depth]=max(FV,[],3);
    
    %% 
    G3=setRange(mean(FV,3));    % mean value in FV along z direction
    G4=setRange(var(FV,0,3));   % variance in FV along z direction
    G5=setRange(maxF);      % max value in FV along z direction
    G6=allinfocus(depth,img);   % AIF image in ImgSeq
    G7=setRange(computeCrossCorr(img,FV)); % corr among img curve & focus curve
    [G8,G9]=computeEigenSV(img);
    [G10,G11]=computeST2D(img);
    [G12,G13]=computeST3D(img);
    
    
end

function output = setRange(In) % scale data between [0 - 1]
    maxIn=max(max(In));
    zeroOneIn=In/maxIn;
    output=round(zeroOneIn*255);
end