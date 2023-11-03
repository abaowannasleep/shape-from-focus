%
%   Mutual-Structure for Joint Filtering
%
%   Parameters:
%   @I0, G0          : Input images, any # of channels
%   @r                  : Window size for filtering, the default is set to 3
%   @eps_I, eps_G       : Control the smoothness of I and G, larger value makes smoother results
%   @lambda_I, lambda_G : Blance the error with repect to I0 and G0 respectively
%   @maxiter            : The maximum iterations for the algorithm
%
%   Note
%   ==========
%   This implementation is the core part of our mutual-structure for joint filtering, so the result
%   of this implementation will be different with that in the corresponding paper because of some
%   pre-processing or post-processing made for the specified task.
%
%   ==========
%   The Code is created based on the method described in the following paper:
%   [1] "Mutual-Structure for Joint Filtering", Xiaoyong Shen, Chao Zhou, Li Xu, Jiaya Jia,
%        IEEE International Conference on Computer Vision (ICCV), 2015
%
%   The code and the algorithm are for non-comercial use only.
%