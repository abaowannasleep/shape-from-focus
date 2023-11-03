function create_monomial_filters_for_structure_tensor_computation3d()
% CREATE_MONOMIAL_FILTERS_FOR_STRUCTURE_TESNOR_COMPUTATION3D Creates 3D monomial filters for structure tensor estimation
%
% create_monomial_filters_for_struncture_tensor_computation3d()
%
% Creates three sets of monomials and stores them in
% 'monomialsForStructureTensor3D.mat'.
%
% INPUT_ARGUMENTS
% N/A
%
% OPTIONAL INPUT ARGUMENTS:
% N/A
%
% OUTPUT ARGUMENTSls

% N/A
%
% See "Advanced filter design" by Knutsson et al for detailed explanation
% of the parameters used.
%
% See also F_weightgensnr, goodsw, bplognorm, krnopt

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

clear coarse intermediate fine
filePath = fileparts(mfilename('fullpath'));
currentPath = pwd;
cd(filePath);
scales = {'coarse','intermediate','fine'};

%%

if ~exist('F_weightgensnr.m','file') || ~exist('goodsw.m','file') || ...
        ~exist('bplognorm.m','file') || ~exist('krnopt.m','file')
    error(['You appear to be missing vital functions for optimizing filters. ',...
        'Please download: https://www.imt.liu.se/edu/courses/TBMI02/code/kerngen.zip'])
end

for k = 1 : length(scales)
    scale = scales{k};
    
    switch scale
        case 'coarse'
            % Filter parameters
            n = 9;
            N = 19;
            spatial_rexp = 2;
            frequency_rexp = -1.0;
            cosexp = 1;
            SNRest = 30;
            DCamp = 1000;
            spacing = [1 1 1];
            fw_amp0 = 30;
            fw_amp1 = 30;
            fw_amp2 = 30;
            
            % Center frequency and bandwidth (in octaves)
            u0 = pi/3;
            B = 2.0;
        case 'intermediate'
            % Filter parameters
            n = 7;
            N = 17;
            spatial_rexp = 2;
            frequency_rexp = -1.0;
            cosexp = 1;
            SNRest = 30;
            DCamp = 1000;
            spacing = [1 1 1];
            fw_amp0 = 30;
            fw_amp1 = 30;
            fw_amp2 = 30;
            
            % Center frequency and bandwidth (in octaves)
            u0 = pi*sqrt(2)/3;
            B = 2.0;
        case 'fine'
            % Filter parameters
            n = 7;
            N = 15;
            spatial_rexp = 2;
            frequency_rexp = -1.0;
            cosexp = 1;
            SNRest = 30;
            DCamp = 1000;
            spacing = [1 1 1];
            fw_amp0 = 30;
            fw_amp1 = 30;
            fw_amp2 = 30;
            
            % Center frequency and bandwidth (in octaves)
            u0 = pi*2/3;
            B = 2.0;
        otherwise
            error('Undefined scale')
    end
    
    fprintf('Optimizing monomials for structure tensor computation in 3d\n');
    
    % Sizes
    spatialSize = [n n n];
    frequencySize = [N N N];
    
    % Frequency weights
    Fw = F_weightgensnr(frequencySize,frequency_rexp,cosexp,SNRest,DCamp);
    
    % Spatial weights
    fw0 = goodsw(spatialSize,spatial_rexp);
    
    % Spatial ideal filter
    fi = wa(spatialSize,0);
    fi = putorigo(fi,1);
    
    % Spatial mask
    fm = wa(spatialSize,0);
    fm = putdata(fm,ones(spatialSize));
    
    % Create coordinate grid
    [Ys, Xs, Zs] = getgrid(wa(frequencySize,1));
    
    Xs = Xs/spacing(1);
    Ys = Ys/spacing(2);
    Zs = Zs/spacing(3);
    
    Rs = sqrt(Xs.*Xs + Ys.*Ys + Zs.*Zs) + eps;
    
    Xh = Xs./Rs;   % h for hat
    Yh = Ys./Rs;
    Zh = Zs./Rs;
    
    % Order 0
    Fi0 = bplognorm(frequencySize, u0, B, spacing);
    
    % Order 1
    Fi1{1} = Fi0.*Xh;
    Fi1{2} = Fi0.*Yh;
    Fi1{3} = Fi0.*Zh;
    
    % Order 2
    Fi2{1}  = Fi0.*Xh.*Xh;
    Fi2{2}  = Fi0.*Xh.*Yh;
    Fi2{3}  = Fi0.*Xh.*Zh;
    Fi2{4}  = Fi0.*Yh.*Yh;
    Fi2{5}  = Fi0.*Yh.*Zh;
    Fi2{6}  = Fi0.*Zh.*Zh;
    
    % Order 0
    fprintf('Order 0------------------------------ \n')
    fw = fw_amp0.*fw0;
    fprintf('Optimizing monomialfilter order 0\n\n')
    [f0, F0] = krnopt(Fi0, Fw, fm, fi, fw);
    
    % Order 1
    fprintf('Order 1------------------------------ \n')
    fw = fw_amp1.*fw0;
    for k = 1: length(Fi1)
        fprintf('Optimizing monomialfilter order 1  %d/%d\n',k,length(Fi1))
        [f1{k}, F1{k}] = krnopt(Fi1{k}, Fw, fm, fi, fw);
    end
    fprintf('\n\n')
    
    % Order 2
    fprintf('Order 2------------------------------ \n')
    fw = fw_amp2.*fw0;
    for k = 1: length(Fi2)
        fprintf('Optimizing monomialfilter order 2 %d/%d\n',k,length(Fi2))
        [f2{k}, F2{k}] = krnopt(Fi2{k}, Fw, fm, fi, fw);
    end
    fprintf('\n\n')
    
    f0 = real(getdata(f0));
    for k = 1: length(Fi1)
        f1{k} = imag(getdata(f1{k}));
    end
    for k = 1: length(Fi2)
        f2{k} = real(getdata(f2{k}));
    end
    
    % Remember filter parameters
    filterInfo.n = n;
    filterInfo.N = N;
    filterInfo.spatial_rexp = spatial_rexp;
    filterInfo.frequency_rexp = frequency_rexp;
    filterInfo.cosexp = cosexp;
    filterInfo.SNRest = SNRest;
    filterInfo.DCamp = DCamp;
    filterInfo.spacing = spacing;
    filterInfo.fw_amp0 = fw_amp0;
    filterInfo.fw_amp1 = fw_amp1;
    filterInfo.fw_amp2 = fw_amp2;
    filterInfo.u0 = u0;
    filterInfo.B = B;
    
    switch scale
        case 'coarse'
            coarse.f0 = f0;
            coarse.f1 = f1;
            coarse.f2 = f2;
            coarse.filterInfo = filterInfo;
        case 'intermediate'
            intermediate.f0 = f0;
            intermediate.f1 = f1;
            intermediate.f2 = f2;
            intermediate.filterInfo = filterInfo;
        case 'fine'
            fine.f0 = f0;
            fine.f1 = f1;
            fine.f2 = f2;
            fine.filterInfo = filterInfo;
    end
end
save monomialsForStructureTensor3D coarse intermediate fine

cd(currentPath)