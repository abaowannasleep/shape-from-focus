%  This is the reference implementation of the fast global image smoother
%  described in the paper:
%
%    Fast Global Image Smoothing based on Weighted Least Squares
%    D. Min, S. Choi, J. Lu, B. Ham, K. Sohn, and M. N. Do,
%    IEEE Trans. Image Processing, vol. no. pp., 2014.

%------------------------------------------------------------------------------

function rDV = WLS2D(Guidance, focusVolume, Lambda, solver_iteration)
    L=length(focusVolume);
    [R,C]=size(focusVolume);
    rDV=zeros(R,C,solver_iteration);
    %rDV=cell(L,1);
    % n=1;
    % m=1;
    % o=1;
    for l=1:L
        %cImg=cImg;
        %cImg=cImg;
        %Guidance=guidance{l};
        %Guidance2=CorGuide{l};
        %     focus_volume=focusVolume{l};
        focus_volume=focusVolume;
        [height, width] = size(focus_volume);
        volume_filtered = zeros(height, width, solver_iteration);
        
        if ~isempty(Guidance) && ~isequal(size(focus_volume), size(Guidance))
            error('Guidance and Data dimentions mismatch');
        end
        
        %     if ismember(l, [2,5,9])
        %         lambda=Lambda(2);
        %     else
        %         lambda=Lambda(1);
        %     end
        lambda=Lambda;
        
        
        %     if l==1
        %         lambda = Lambda(2);   % 50
        %     elseif l==2
        %            lambda = Lambda(4);   % 10
        %     elseif l==3
        %            lambda = Lambda(5);    % 5
        %     elseif l==4
        %            lambda = Lambda(1);    % 60
        %     elseif l==5
        %            lambda = Lambda(3);   % 15
        %     elseif l==6
        %            lambda = Lambda(3);   % 15
        %     elseif ismember(l, [7,8])
        %            lambda = Lambda(4);   % 10
        %     elseif l==9
        %            lambda = Lambda(6);   % 2
        %     end
        
        
        MaxDif2=findMaxDif2(Guidance);
        BLFKernel = prepareBLFKernel(MaxDif2);
        
        
        % BLFKernelI = prepareBLFKernel(sigma);
        
        zeros_width = zeros(width, 1);
        zeros_height = zeros(height, 1);
        
        if isempty(Guidance)
            Guidance = focus_volume;
            %jointFiltering = 0;
            %else
            %jointFiltering = 1;
        end
        
        lambda_in = lambda;
        
        for iter=1:solver_iteration
            
            
            %lambda_in = 1.5 * lambda * 4 ^ (solver_iteration - iter) / (4 ^ solver_iteration - 1);
            
            % -- Filter levels individually
            % for each row
            % for each column
            for j=1:width
                a_vec = zeros_height;
                c_vec = zeros_height;
                for i=2:height
                    % compute bilateral weight for all channels
                    %color_diff = (focus_volume(i,j,k) - focus_volume(i-1,j,k))^2;
                    %color_diff = abs(Guidance(i,j,k) - Guidance(i-1,j,k));
                    %a_vec(i) = -lambda_in * BLFKernelI(round(color_diff)+1);
                    var_diff = abs(Guidance(i,j) - Guidance(i-1,j));
                    %var_dsum = ((abs(Guidance(i,j,k) - Guidance(i-1,j,k))+1)^2)/(Guidance(i,j,k) + Guidance(i-1,j,k)+1);
                    %var_sum = Guidance(i,j,k) + Guidance(i-1,j,k);
                    %a_vec(i) = -lambda_in * BLFKernelI(round(var_dsum+1));
                    a_vec(i) = -lambda_in * BLFKernel(round(var_diff+1));
                    %a_vec(i) = -lambda_in * BLFKernelI(round(1/(var_diff+1)+1));
                    %                        misco2(m)=round(1/(var_sum+1)+1);
                    %                        m=m+1;
                    %a_vec(i) = -lambda_in * exp(-3*Guidance2(i,j,2)/2);
                end
                c_vec(1:end-1) = a_vec(2:end);
                b_vec = 1 - a_vec - c_vec;
                focus_volume(:,j) = solve_tridiagonal(focus_volume(:,j)', height, a_vec, b_vec, c_vec);
            end
            
            for i=1:height
                a_vec = zeros_width;
                c_vec = zeros_width;
                for j=2:width
                    % compute bilateral weight for all channels
                    
                    var_diff = abs(Guidance(i,j) - Guidance(i,j-1));
                    
                    a_vec(j) = -lambda_in * BLFKernel(round(var_diff+1));
                    
                end
                c_vec(1:end-1) = a_vec(2:end);
                b_vec = 1 - a_vec - c_vec;
                focus_volume(i,:) = solve_tridiagonal(focus_volume(i,:)', width, a_vec, b_vec, c_vec);
            end
            %--------------------------------------------------------------
            lambda_in = lambda_in/2;
            
            volume_filtered(:,:,iter) = focus_volume;
        end
        rDV= volume_filtered;
    end
    
end

%------------------------------------------------------------------------------------------
% function BLFKernelI = prepareBLFKernel(sigma)
% %PREPAREBLFKERNEL Summary of this function goes here
% %   Detailed explanation goes here
% 	MaxSizeOfFilterI = 111165025;
%    % MaxSizeOfFilterI = MaxDif+2;       function BLFKernelI = prepareBLFKernel(sigma, MaxDif)
%    BLFKernelI= exp(-sqrt(0:(MaxSizeOfFilterI-1)) / sigma);
%    %BLFKernelI = exp(-(0:(MaxSizeOfFilterI-1)) / sigma);
% end
%------------------------------------------------------------------------------------------
function BLFKernelI = prepareBLFKernel(MaxSizeOfFilterI)
    %PREPAREBLFKERNEL Summary of this function goes here
    %   Detailed explanation goes here
    MaxSizeOfFilter=round(MaxSizeOfFilterI+1);
    x=0:MaxSizeOfFilter;
    BLFKernelI = exp(-3*x/MaxSizeOfFilter);
end

%--------------------------------------------------------------------------------------------
function [x] = solve_tridiagonal(x_in, N, a, b, c)
    x = x_in;
    c(1) = c(1) / b(1);
    x(1) = x(1) / b(1);
    
    % loop from 1 to N-1 inclusive
    for n=2:N
        m = 1 / (b(n) - a(n) * c(n-1));
        c(n) = c(n) * m;
        x(n) = (x(n) - a(n) * x(n-1)) * m;
    end
    
    
    % loop from N-2 to 0 inclusive
    for n=N-1:-1:1
        x(n) = x(n) - c(n) * x(n+1);
    end
    
end

%
function MaxDif=findMaxDif2(A)   % max neighborhood diff of 2d image considering all 4 directions
    P=padarray(A,[1,1],'replicate');
    right=P(2:end-1,3:end);
    left=P(2:end-1,1:end-2);
    above=P(1:end-2,2:end-1);
    below=P(3:end,2:end-1);
    DifL=A-left;
    DifR=A-right;
    DifA=A-above;
    DifBe=A-below;
    maxL=max(max(max(abs(DifL))));
    maxR=max(max(max(abs(DifR))));
    maxA=max(max(max(abs(DifA))));
    maxBe=max(max(max(abs(DifBe))));
    MaxDif=max([maxA,maxBe,maxL,maxR]);
    
end