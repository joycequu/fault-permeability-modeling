% this was used to generate the plots

folderPath = '/Users/joycequ/Documents/UROP/urop24-joyce/predict/tornado_plots/1layer';

% Specify the filename (including the .mat extension)
filename = '1layer_perm_mean_output.mat';

% Construct the full file path
fullFilePath = fullfile(folderPath, filename);

% Load the .mat file
data = load(fullFilePath);

% data = load('faults_perm_mean_output.mat', '-mat');
output_data = data.faults_perm_mean;

% variables = transpose(output_data(2:11, 1));
kxx_effect = zeros(1, 10);
kyy_effect = zeros(1, 10);
kzz_effect = zeros(1, 10);

variables = {'vcl_fw', 'vcl_hw', 'fdip', 'zf', 'zmax'};

kxx_effect_lower = zeros(1, 5);
kxx_effect_upper = zeros(1, 5);

kyy_effect_lower = zeros(1, 5);
kyy_effect_upper = zeros(1, 5);

kzz_effect_lower = zeros(1, 5);
kzz_effect_upper = zeros(1, 5);

baseline = cell2mat(output_data(1, 2));

for i = 1 : 10
    var_effect = cell2mat(output_data(i + 1, 2));
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

kxx_effect_lower = log10(kxx_effect_lower./(milli*darcy));
kxx_effect_upper = log10(kxx_effect_upper./(milli*darcy));

kyy_effect_lower = log10(kyy_effect_lower./(milli*darcy));
kyy_effect_upper = log10(kyy_effect_upper./(milli*darcy));

kzz_effect_lower = log10(kzz_effect_lower./(milli*darcy));
kzz_effect_upper = log10(kzz_effect_upper./(milli*darcy));

% Sort the values based on the lower change
% Sort the higher values and the names arrays using the same indices
% [kyy_effect_lower, ind] = sort(kyy_effect_lower,'descend');
% kyy_effect_upper = kyy_effect_upper(ind);
% variables_sorted = variables(ind);

% Create a figure and plot the low and high horizontally
figure;
h1 = barh(kxx_effect_lower, 'b');
hold on;
h2 = barh(kxx_effect_upper,'r');
bh1 = get(h1,'BaseLine');
set(bh1,'BaseValue',baseline(2));
title('Sensitivity Analysis (kxx)');
% yticks(1:numel(variables_sorted));
% yticklabels(variables_sorted);
yticks(1:numel(variables));
yticklabels(variables);
xlabel('Mean Permeability');
legend([h1 h2],'Negative Perturbation','Positive Perturbation')

text(h1.XEndPoints, h1.YEndPoints, string(h1.YEndPoints), 'Horiz','right', 'Vert','middle')
text(h2.XEndPoints, h2.YEndPoints, string(h2.YEndPoints), 'Horiz','left', 'Vert','middle')

saveas(gcf, 'local_sens_tornado_sorted.png');