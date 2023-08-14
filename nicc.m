function niccErr = nicc(fm,fp)
fp = double(abs(fp - fp(:,:,1)));
std1 = std(fm, 0, 3);
std2 = std(fp, 0, 3);
m1 = mean(fm, 3);
m2 = mean(fp, 3);
[~, ~, p] = size(fp);
temp = (fm-m1).*(fp-m2);
cov_mp = sum(temp, 3)/(p-1);
niccErr = cov_mp ./ (std1 .* std2);
end