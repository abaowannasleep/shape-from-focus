function [t11, t12, t13, t22, t23, t33] = compute_structure_tensor3d(input, varargin)
% COMPUTE_STRUCTURE_TENSOR3D Computes the structure of the provided image
%
% function [t11 t12 t13 t22 t23 t33] =
% compute_structure_tensor3d(input, average, normalize)
%
% INPUT ARGUMENTS
% input                     - Input data to estimate the structure tensor for
%
% OPTIONAL INPUT ARGUMENTS
% 'scale'                   - Scale of filters to apply, coarse, intermediate or
%                             fine (intermediate)
% 'average'                 - Set to true to average the estimated tensor (false)
% 'sizeAveragingFilter'     - Spatial size of averaging filter (7)
% 'sigmaAveragingFilter'    - Sigma of averaging filter, defined in the FD
%                             domain (1.0)
% 'normalize'               - Set to true to normalize the estimated tensor (false)
% 'mode'                    - 'quadrature'/'monomials'
%
% OUTPUT ARGUMENTS
% t11                       - Tensor element 1,1
% t12                       - Tensor element 1,2
% t13                       - Tensor element 1,3
% t22                       - Tensor element 2,2
% t23                       - Tensor element 2,3
% t33                       - Tensor element 3,3
%
% Please see "Representing local structure using tensors" by H Knutsson
% and "Representing local structure using tensors II" by H Knutsson et al.
% for theory behind structure tensors and the use of quadrature and/or
% monomials for computing them. The theory is also discussed extensively in
% "Signal processing for computer vision" by G Granlund and H Knutsson.

% Copyright (c) 2012 Daniel Forsberg
% danne.forsberg@outlook.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Set up parameter default values
scale = 'intermediate';
average = false;
sizeAveragingFilter = 5;
sigmaAveragingFilter = 1.0;
normalize = false;
mode = 'quadrature';

%% Overwrites default parameter
for k=1:2:length(varargin)
    eval([varargin{k},'=varargin{',int2str(k+1),'};']);
end

% Initialize tensor elements
t11 = zeros(size(input));
t12 = zeros(size(input));
t13 = zeros(size(input));
t22 = zeros(size(input));
t23 = zeros(size(input));
t33 = zeros(size(input));

if strcmp(mode,'quadrature')
    % Load filters
    if ~exist('quadratureFiltersForStructureTensor3D.mat','file')
        disp('Quadrature filters needed for compute_structure_tensor3d are not available.')
        answer = input('Would you like to download these filters? [y]/[n] ','s');
        if strcmpi(answer,'y') || strcmpi(answer,'yes')
            folder = fileparts(mfilename('fullpath'));
            currentFolder = pwd;
            cd(folder)
            urlwrite(...
                'https://github.com/fordanic/tensor-processing/blob/master/quadratureFiltersForStructureTensor3D.mat',...
                'quadratureFiltersForStructureTensor3D.mat')
            cd(currentFolder)
        end
    end
    load quadratureFiltersForStructureTensor3D
    
    % Select filters
    eval(['qFilt = ',scale,'.f;']);
    eval(['m11 = ',scale,'.m11;']);
    eval(['m12 = ',scale,'.m12;']);
    eval(['m13 = ',scale,'.m13;']);
    eval(['m22 = ',scale,'.m22;']);
    eval(['m23 = ',scale,'.m23;']);
    eval(['m33 = ',scale,'.m33;']);
    
    for k = 1 : 6
        q = imfilter(input,qFilt{k},'conv','same','replicate');
        
        % Estimate T
        Aq = abs(q);
        t11 = t11 + Aq*m11{k};
        t12 = t12 + Aq*m12{k};
        t13 = t13 + Aq*m13{k};
        t22 = t22 + Aq*m22{k};
        t23 = t23 + Aq*m23{k};
        t33 = t33 + Aq*m33{k};
    end
elseif strcmp(mode,'monomials')
    % Load filters
    if ~exist('monomialsForStructureTensor3D.mat','file')
        disp('Monomial filters needed for compute_structure_tensor3d are not available.')
        answer = input('Would you like to download these filters? [y]/[n] ','s');
        if strcmpi(answer,'y') || strcmpi(answer,'yes')
            folder = fileparts(mfilename('fullpath'));
            currentFolder = pwd;
            cd(folder)
            urlwrite(...
                'https://github.com/fordanic/tensor-processing/blob/master/monomialsForStructureTensor3D.mat',...
                'monomialsForStructureTensor3D.mat')
            cd(currentFolder)
        end
    end
    load monomialsForStructureTensor3D
    
    % Select filters
    eval(['mFilt1 = ',scale,'.f1;']);
    eval(['mFilt2 = ',scale,'.f2;']);
    
    for k = 1 : length(f1)
        q1{k} = imfilter(input,mFilt1{k},'conv','same','replicate');
    end
    for k = 1 : length(f2)
        q2{k} = imfilter(input,mFilt2{k},'conv','same','replicate');
    end
    
    t11  = q1{1}.*q1{1};
    t12  = q1{1}.*q1{2};
    t13  = q1{1}.*q1{3};
    t22  = q1{2}.*q1{2};
    t23  = q1{2}.*q1{3};
    t33  = q1{3}.*q1{3};
    
    t11 = t11 + ...
        q2{1}.*q2{1} + q2{2}.*q2{2} + q2{3}.*q2{3};
    
    t12 = t12 + ...
        q2{1}.*q2{2} + q2{2}.*q2{4} + q2{3}.*q2{5};
    
    t13 = t13 + ...
        q2{1}.*q2{3} + q2{2}.*q2{5} + q2{3}.*q2{6};
    
    t22 = t22 + ...
        q2{2}.*q2{2} + q2{4}.*q2{4} + q2{5}.*q2{5};
    
    t23 = t23 + ...
        q2{2}.*q2{3} + q2{4}.*q2{5} + q2{5}.*q2{6};
    
    t33 = t33 + ...
        q2{3}.*q2{3} + q2{5}.*q2{5} + q2{6}.*q2{6};
    
    Tnorm = (t11.^2 + 2*t12.^2 + 2*t13.^2 + t22.^2 + 2*t23.^2 + t33.^2).^(1/4) + eps;
    
    t11 = t11./Tnorm;
    t12 = t12./Tnorm;
    t13 = t13./Tnorm;
    t22 = t22./Tnorm;
    t23 = t23./Tnorm;
    t33 = t33./Tnorm;
else
    error('Unknown computation mode')
end
if average
    Tcert = sqrt(t11.^2 + 2*t12.^2 + 2*t13.^2 + t22.^2 + 2*t23.^2 + t33.^2);
    
    t11 = averaging3d(t11, Tcert, sizeAveragingFilter, sigmaAveragingFilter);
    t12 = averaging3d(t12, Tcert, sizeAveragingFilter, sigmaAveragingFilter);
    t13 = averaging3d(t13, Tcert, sizeAveragingFilter, sigmaAveragingFilter);
    t22 = averaging3d(t22, Tcert, sizeAveragingFilter, sigmaAveragingFilter);
    t23 = averaging3d(t23, Tcert, sizeAveragingFilter, sigmaAveragingFilter);
    t33 = averaging3d(t33, Tcert, sizeAveragingFilter, sigmaAveragingFilter);
end

if normalize
    maxTNorm = max(vec(sqrt(t11.^2 + 2*t12.^2 + 2*t13.^2 + t22.^2 + 2*t23.^2 + t33.^2 + eps)));
    
    t11 = t11/maxTNorm;
    t12 = t12/maxTNorm;
    t13 = t13/maxTNorm;
    t22 = t22/maxTNorm;
    t23 = t23/maxTNorm;
    t33 = t33/maxTNorm;
end
