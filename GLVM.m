function FMs = GLVM(images, WSize, smoothSize)
%
% Z. M. He
% function GLVM : GLVM focus measure.
%
% Input:
%     images: n x m x p matrix; p: number of images; 
%     Wsize: window size for focus measure;
%     smoothSize: smooth window size;
% Output:
%     FM: measured focus;
%
[height,width,imgNum]=size(images);
FMs = zeros(height,width,imgNum);
MEANF = fspecial('average',[WSize WSize]);
for k = 1:imgNum
    Image = images(:,:,k);
    U = imfilter(Image, MEANF, 'replicate');
    FM= (Image-U).^2;
%     FM = imfilter(FM, MEANF, 'replicate'); 
    FMs(:, :, k) = FM;
end
if smoothSize
    FMs = smooth3(FMs,'box',smoothSize);
end
end