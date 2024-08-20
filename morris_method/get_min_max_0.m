data = load('ee_result_500.mat', '-mat');
ee = data.ee;

mean_vals = zeros(4, 6, 3);
std_vals = zeros(4, 6, 3);

for d = 1 : 3
    mean_vals(:, :, d) = ee(:, :, d*2-1);
    std_vals(:, :, d) = ee(:, :, d*2);
end

min_mean = min(mean_vals, [], 'all');
max_mean = max(mean_vals, [], 'all');
min_std = min(std_vals, [], 'all');
max_std = max(std_vals, [], 'all');
