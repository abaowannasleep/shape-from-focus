function [me_err] = me(mat1,mat2,mask)
%SSE 此处显示有关此函数的摘要
%   此处显示详细说明
err = abs((mat1-mat2).*mask);
me_err = max(max(err));
end

