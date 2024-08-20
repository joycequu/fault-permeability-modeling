% this was not used to generate the plots

data = load('faults_perm_mean_output.mat', '-mat');
output_data = data.faults_perm_mean;

% variables = transpose(output_data(2:11, 1));
kxx_effect = zeros(1, 10);
kyy_effect = zeros(1, 10);
kzz_effect = zeros(1, 10);

variables = {'vcl_fw', 'vcl_hw', 'sfdip', 'szf', 'szmax'};
kxx_effect_lower = zeros(1, 5);
kxx_effect_upper = zeros(1, 5);
kyy_effect_lower = zeros(1, 5);
kyy_effect_upper = zeros(1, 5);
kzz_effect_lower = zeros(1, 5);
kzz_effect_upper = zeros(1, 5);

baseline = cell2mat(output_data(1, 2));

% for i = 1 : 10
%     var_effect = cell2mat(output_data(i + 1, 2)) - baseline;
%     kxx_effect(i) = var_effect(1);
%     kyy_effect(i) = var_effect(2);
%     kzz_effect(i) = var_effect(3);
% end

for i = 1 : 10
    var_effect = cell2mat(output_data(i + 1, 2)) - baseline;
    if mod(i, 2) == 1
        kxx_effect_lower((i+1)/2) = var_effect(1);
        kyy_effect_lower((i+1)/2) = var_effect(2);
        kzz_effect_lower((i+1)/2) = var_effect(3);
    else
        kxx_effect_upper(i/2) = var_effect(1);
        kyy_effect_upper(i/2) = var_effect(2);
        kzz_effect_upper(i/2) = var_effect(3);
    end
end

% Plot tornado plot
figure;
barh(kyy_effect_lower, 'b');
hold on;
barh(kyy_effect_upper, 'r');
yticks(1:numel(variables));
yticklabels(variables);

xlabel('Effect on Mean Permeability (kyy)');

title('Local Sensitivity Analysis Plot (Small)');
% legend('Location', 'Best');
grid on;