function setup_tensor_processing_repository()
% SETUP_TENSOR_PROCESSING_REPOSITORY Call this function to setup this repository
%
% setup_tensor_processing_repository
%
% INPUT ARGUEMNTS
% N/A
%
% OPTIONAL INPUT ARGUEMNTS
% N/A
%
% OUTPUT ARGUEMNTS
% N/A

% Copyright (c) 2015 Daniel Forsberg
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

% Get folder of current file and go to folder
folder = fileparts(mfilename('fullpath'));
currentFolder = pwd;
cd(folder)

%% Download external dependencies if needed

% Load external dependencies and check which are missing
load externalDependenciesTensorProcessing
missingExternalDependencies = [];
for k = 1 : length(fileNames)
    if ~exist(fileNames{k},'file')
        missingExternalDependencies(end + 1) = k;
    end
end

% Ask missing files should be downloaded
if ~isempty(missingExternalDependencies)
    disp('The following files, which are necessary for the tensor processing repository are missing:')
    for k = 1 : length(missingExternalDependencies)
        disp(fileNames{k})
    end
    answer = input('Would you like to download them? [y]/[n] ','s');
    
    % If yes, download missing files
    if strcmp(lower(answer),'y') || strcmp(lower(answer),'yes')
        if ~exist([folder,filesep,'external-dependencies'],'dir')
            mkdir('external-dependencies')
        end
        cd('external-dependencies')
        for k = 1 : length(missingExternalDependencies)
            disp(['Downloading: ',fileNames{missingExternalDependencies(k)}])
            urlwrite(urls{missingExternalDependencies(k)},...
                fileNames{missingExternalDependencies(k)});
        end
        addpath(pwd)
    end
end

%% Add external dependencies path
if exist([folder,filesep,'external-dependencies'],'dir')
    addpath([folder,filesep,'external-dependencies'])
end

% Go back to original folder
cd(currentFolder)

