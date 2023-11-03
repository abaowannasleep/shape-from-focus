function newDepth = computeGuidedDepthBalls(d,g)
    %% guidance image
    d1 = g;
    
    %% Initial depth
    d2 = d;
    
    %% Digital Photography with Flash and No-Flash Image Pairs
    d3=bfilter2(d,g);
    
    %% An application of MRF to range sensing
    d4=mrf(d,g,8,0.1);   
    
    %% Joint bilateral upsampling
    d5 = JBU(d,g,1,5,0.1,5);
    
    %% A noise-aware filter for real-time depth upsampling
    d6 = NAFilter(d,g,1,5,10,2,3,0.5,10);
    
    %% Image and sparse laser fusion for dense scene reconstruction
    d7 = LSF(d,g,5,5);
    
    %% Depth video enhancement based on weighted mode filtering
    d8 = WMoF(d,g,3,10,10,3);
    
    %% Guided Image Filtering
    d9=imguidedfilter(d,g);
    
    %% Guided depth enhancement via anisotropic diffusion
    d10 = AD(d,g,1,0.1);
    
    %% Constant time weighted median filtering for stereo matching and beyond
    eps = 0.01^2;
    r = ceil(max(size(g,1), size(g,2)) / 40);
    d11 = weighted_median_filter(d,repmat(g,[1 1 3]), 0:255, r, eps); % r=5 ?
    
    %% 100+ times faster weighted median filter
    d12=jointWMF(d,g,10,40,256,256,1,'exp');
    
    %% Weighted Guided Image Filtering
    radius=4;
    d13=WeightedGuidedImageFilter(d, g, radius, 0.5, 1);
%     %% Rolling Guidance Filter
%     d6 = RollingGuidanceFilter(d,10,0.2,10);

%% Gradient Domain Guided Image Filtering
    radius=4;
    lambda=0.01; % 10
%     def=d./max(max(d)); gef=g./max(max(g));
    d14 = gguidedfilter(g, d, radius, lambda);
%     d8=d8.*max(max(d));
    %     D01=D-min(min(D));  % D is ground truth depth needed for d11 scaling
    %     d11_01=d11-min(min(d11));
    %     d11=(max(max(D01))/max(max(d11_01))).*d11_01+min(min(D));
    
    %% Mutual Structure for Joint Filtering
    eps_I = 5e-4;
    eps_G = 10e-4;
    lambda_I = 100;    % for target image
    lambda_G = 100;    % for guidance image
    maxiter = 20;
    %     r = 2; % it was bad
    r = 4;
    [~,d15]= mutual_structure_joint_filtering(repmat(g,[1 1 3]), medfilt2(d,[3 3]), r, eps_I, eps_G, lambda_I, lambda_G, maxiter);  
    
    %% Fast domain decomposition for global image smoothing (WL1)
    sigma_r=0.02;
    lambda=50;
    d16 = fwl1(repmat(d,[1 1 3]), repmat(g,[1 1 3]), lambda, 5,5, sigma_r);
    d16 = reshape(d16,size(repmat(d,[1 1 3]),1),size(repmat(d,[1 1 3]),2),size(repmat(d,[1 1 3]),3));
    d16 = d16(:,:,2);
    
    %% Fast domain decomposition for global image smoothing (WLS)
    sigma_r=0.2;
    lambda=3^2;
    d17 = fwls(repmat(g,[1 1 3]), repmat(d,[1 1 3]), sigma_r,lambda, 5, 1.2);
    d17 = reshape(d17,size(repmat(d,[1 1 3]),1),size(repmat(d,[1 1 3]),2),size(repmat(d,[1 1 3]),3));
    d17 = d17(:,:,2);
    
    %% Robust Guided Image Filtering using Nonconvex Potentials
    [m, n] = size(d);
    u0=ones(m,n);
    nei= 0; % better than 1
    lambda = 20^2;           % for real
    mu = 500;              % bw for static guidance
    nu = 200;              % bw for dynamic guidance
    step=10;
    issparse = true;
    d18 = sdfilter(g,u0,d,nei,lambda,mu,nu,step,issparse);
    
    %% Effective Guided Image Filtering for Contrast Enhancement
    r = 6;
    eps = 500;  
%     def=d./max(max(d)); gef=g./max(max(g));
    [q,beta]=egif(g, d, r, eps);
    d19=(g - q).* beta  + q;
%     d13=d13.*max(max(d));
    
    %% Mutually Guided Image Filtering
    %     alpha_t = .001;% can be tuned
    %     alpha_r = 0;
    %     N_iter =10;
    %     mode = 2;    % mode = 2; % dynamic only
    %     [d14,~] = muGIF(d,g,alpha_t,alpha_r,N_iter,mode);
    %
    alpha_t = 1;% can be tuned
    alpha_r = 0;
    N_iter =10;
    mode = 1;   % mode = 1; % static/dynamic
    [d20,~] = muGIF(d,g,alpha_t,alpha_r,N_iter,mode);
    %
    alpha_t = .005;% can be tuned
    alpha_r = 0.02;
    N_iter =10;
    mode = 0;   %mode = 0 dynamic/dynamic
    [d21,~] = muGIF(d,g,alpha_t,alpha_r,N_iter,mode);
    
    %%
    ofst = 5;
    Depth=cell(21,1);
    newDepth=cell(21,1);
    for i=1:21
        Depth{i,1}=eval(sprintf('d%d',i));
        newDepth{i,1}=Depth{i,1}(ofst+1:end-ofst,ofst+1:end-ofst);
    end
end