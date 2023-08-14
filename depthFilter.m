function mask = depthFilter(z)
kernel = [[1, 1, 1];[1, -8, 1];[1, 1, 1]];
mask = zeros(size(z));
z_filter = 20*abs(imfilter(z, kernel, 'replicate'));
mask(z_filter>6) = 1;
% hole fill
% se = strel('square', 3);
% bmask = imdilate(mask, se);
bmask = imfill(mask, 'holes');
if sum(sum(bmask)) > 3*sum(sum(mask))
    mask = 1-mask;
else
    mask = 1-bmask;
end
% bmask = imfill(mask, 'holes');
% mask = 1-imfill(mask, 'holes');
end