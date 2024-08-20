% this was used to generate the plots

data = load('sc_c_2layer_perm_mean_output.mat', '-mat');
output_data = data.faults_perm_mean;

% variables = transpose(output_data(2:11, 1));
kxx_effect = zeros(1, 12);
kyy_effect = zeros(1, 12);
kzz_effect = zeros(1, 12);

variables = {'vcl (fw_1)', 'vcl (fw_2)', 'vcl (hw)', 'fdip', 'zf', 'zmax'};

kxx_effect_lower = zeros(1, 6);
kxx_effect_upper = zeros(1, 6);

kyy_effect_lower = zeros(1, 6);
kyy_effect_upper = zeros(1, 6);

kzz_effect_lower = zeros(1, 6);
kzz_effect_upper = zeros(1, 6);

kxx_baseline = log10(output_data(1, 1)/(milli*darcy));
kyy_baseline = log10(output_data(1, 2)/(milli*darcy));
kzz_baseline = log10(output_data(1, 3)/(milli*darcy));

for i = 1 : 12
    if mod(i, 2) == 1
        kxx_effect_lower((i+1)/2) = output_data(i+1, 1);
        kyy_effect_lower((i+1)/2) = output_data(i+1, 2);
        kzz_effect_lower((i+1)/2) = output_data(i+1, 3);
    else
        kxx_effect_upper((i/2))= output_data(i+1, 1);
        kyy_effect_upper((i/2)) = output_data(i+1, 2);
        kzz_effect_upper((i/2)) = output_data(i+1, 3);
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

% for kxx
% h1 = barh(kxx_effect_lower, 'b');
% hold on;
% h2 = barh(kxx_effect_upper,'r');
% bh1 = get(h1,'BaseLine');
% set(bh1,'BaseValue', kxx_baseline);
% title('Sensitivity Analysis (2 layers, sc/c, kxx)');

% for kyy
% h1 = barh(kyy_effect_lower, 'b');
% hold on;
% h2 = barh(kyy_effect_upper,'r');
% bh1 = get(h1,'BaseLine');
% set(bh1,'BaseValue', kyy_baseline);
% title('Sensitivity Analysis (2 layers, sc/c, kyy)');

% % for kzz
h1 = barh(kzz_effect_lower, 'b');
hold on;
h2 = barh(kzz_effect_upper,'r');
bh1 = get(h1,'BaseLine');
set(bh1,'BaseValue', kzz_baseline);
title('Sensitivity Analysis (2 layers, sc/c, kzz)');

% yticks(1:numel(variables_sorted));
% yticklabels(variables_sorted);
yticks(1:numel(variables));
yticklabels(variables);
xlabel('Mean Permeability');
legend([h1 h2],'Negative Perturbation','Positive Perturbation')

text(h1.XEndPoints, h1.YEndPoints, string(h1.YEndPoints), 'Horiz','right', 'Vert','middle')
text(h2.XEndPoints, h2.YEndPoints, string(h2.YEndPoints), 'Horiz','left', 'Vert','middle')

% caption = 'FW: sand, clay; HW: clay';
% text(0.5, 0.02, caption, 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'blue', 'FontWeight', 'bold');

saveas(gcf, 'sc_c_kzz_tornado.png');