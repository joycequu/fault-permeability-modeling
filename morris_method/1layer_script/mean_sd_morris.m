% 1 - vcl_fw
% 2 - vcl_hw
% 3 - sfdip
% 4 - szf
% 5 - szmax
% (referred to as v below)
 
data = load('all_ee_twenty.mat', '-mat');
output = data.all_ee;

% output (3d array)
% 1000 for 1000 runs
% 5 for the parameters
% 3 for kxx, kyy, kzz
 
ee = zeros(5, 6);
% 5 for 5 parameters
% 6: mean, variance fore each kxx/kyy
 
for d = 1 : 3
    % d=1 kxx; d=2 kyy; d=3 kzz
    for v = 1 : 5
        k_numeric = output(:, v, d);
        k_abs = abs(k_numeric);
        miu = mean(k_abs);
        ee(v, 2*d-1) = log10(miu./(milli*darcy));
        var = 0;
        % calculate variance
        for i = 1 : 1000
            var = var + (k_numeric(i) - miu)^2;
        end
        var = var / (20 - 1);
        ee(v, 2*d) = log10(var./(milli*darcy)^2);
    end
end
 
save('ee_result.mat', 'ee');
