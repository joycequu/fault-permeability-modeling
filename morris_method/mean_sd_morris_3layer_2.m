% 1 - vcl_fw
% 2 - vcl_hw
% 3 - sfdip
% 4 - szf
% 5 - szmax
% (referred to as v below)
 
data = load('scaled_fault_perm_3layer.mat', '-mat');
output = data.all_ee;

% output (3d array)
% 1000 for 1000 runs
% 5 for the parameters
% 3 for kxx, kyy, kzz

baseline_selection = 500;
% (500, 4, 7, 3)
% 500 for num of baseline selections
% 4 for the sand/clay combos
% 7 for baseline + 6 parameters
% 3 for kxx, kyy, kzz

ee = zeros(4, 6, 6);
% c = 1, sand (fw1), clay (fw2), clay (hw)
% c = 2, clay (fw1), sand (fw2), clay (hw)
% c = 3, sand (fw1), clay (fw2), sand (hw)
% c = 4, clay (fw1), sand (fw2), sand (hw)

% 6 for 6 parameters
% 6: mean, variance fore each kxx/kyy

% THIS IS FOR WHEN INPUT IS 500 * 4 * 7 * 3
for c = 1 : 4
    for d = 1 : 3
        % d=1 kxx; d=2 kyy; d=3 kzz
        for v = 2 : 7
            % v for each parameter (v = 1, baseline)
            k_numeric = output(:, c, v, d) - output(:, c, 1, d);
            miu = mean(abs(k_numeric));
            % need to take abs?
            ee(c, v-1, 2*d-1) = miu;
            var = 0;
            % calculate variance
            for i = 1 : baseline_selection
                var = var + (k_numeric(i) - miu)^2;
            end
            ee(c, v-1, 2*d) = sqrt(var / (baseline_selection - 1));
        end
    end
end
 
save('scaled_ee_result.mat', 'ee');
