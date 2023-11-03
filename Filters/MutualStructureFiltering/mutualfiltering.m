function G = mutualfiltering(depth,rgb)
% rgb = im2double(imread('testdata/rgb_depth_24/rgb_input.png'));
% depth = im2double(imread('testdata/rgb_depth_24/depth_input.png'));

eps_I = 1e-4;
eps_G = 5e-4;
r = 2;
lambda_I = 1;
lambda_G = 30;
maxiter = 15;

[I, G] = mutual_structure_joint_filtering(rgb, medfilt2(depth,[3 3]), r, eps_I, eps_G, lambda_I, lambda_G, maxiter);
G=double(G);

% figure,imshow(I)
% figure,imshow(G)
end