% this was used to generate the plots

% clear everything and previously generate plots
clear
close all force

data = load('2layers_1k_faults_perm_mean_output.mat', '-mat');
total_data = data.faults_perm_mean;

total_names = {'sc_c', 'sc_s', 'cs_c', 'cs_s'};
graph_names = {'SC/C', 'SC/S', 'CS/C', 'CS/S'};

cmap = copper(12);

for g = 1 : 4
    output_data = total_data(13*(g-1)+1 : 13*g, :);
    output_name = total_names{g};
    
    kxx_effect = zeros(1, 12);
    kyy_effect = zeros(1, 12);
    kzz_effect = zeros(1, 12);

    var_names = {'V_{cl} (FW_1)', 'V_{cl} (FW_2)', 'V_{cl} (HW)', 'f_{dip}', 'z_f', 'z_{max}'};

    kxx_effect_lower = zeros(1, 6);
    kxx_effect_upper = zeros(1, 6);

    kyy_effect_lower = zeros(1, 6);
    kyy_effect_upper = zeros(1, 6);

    kzz_effect_lower = zeros(1, 6);
    kzz_effect_upper = zeros(1, 6);

    % change: already took log, so just avg here
    % kxx_baseline = log10(output_data(1, 1)/(milli*darcy));
    % kyy_baseline = log10(output_data(1, 2)/(milli*darcy));
    % kzz_baseline = log10(output_data(1, 3)/(milli*darcy));
    kxx_baseline = output_data(1,1);
    kyy_baseline = output_data(1,2);
    kzz_baseline = output_data(1,3);

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

    % Sort the values based on the lower change
    % Sort the higher values and the names arrays using the same indices
    % [kyy_effect_lower, ind] = sort(kyy_effect_lower,'descend');
    % kyy_effect_upper = kyy_effect_upper(ind);
    % variables_sorted = variables(ind);
    
    % Create a figure and plot the low and high horizontally
    % Plot for kxx-direction
    kxx_fig = figure;
    h1 = barh(kxx_effect_lower, 'FaceColor', cmap(2, :));
    hold on;
    h2 = barh(kxx_effect_upper, 'FaceColor', cmap(10, :));
    bh1 = get(h1,'BaseLine');
    set(bh1,'BaseValue', kxx_baseline);
    graph_title = ['Sensitivity Analysis (2 layers, 1k simulations)' ' for ' graph_names{g} ' in kxx-direction'];
    title(graph_title, 'FontSize', 13);


    % simply labeling
    yticks(1:numel(var_names));
    yticklabels(var_names);
    ax = gca;
    ax.YAxis.FontSize = 14; % for readability purpose
    
    % absolute scaling
    xlim([-7 2.5]);
    xticks([-7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5]);
    xlabel('Mean Permeability');
    ax = gca;
    ax.XAxis.FontSize = 12; % for readability purpose
    ax.XLabel.FontSize = 12;
    legend([h1 h2],'Negative Perturbation','Positive Perturbation')
    % text(h1.XEndPoints, h1.YEndPoints, string(h1.YEndPoints), 'Horiz','right', 'Vert','middle')
    % text(h2.XEndPoints, h2.YEndPoints, string(h2.YEndPoints), 'Horiz','left', 'Vert','middle')
    
    kxx_fig.Position = [100, 100, 800, 400];

    filename = ['abs_' output_name '_kxx_2layer_1k.png'];
    exportgraphics(kxx_fig, filename, 'Resolution', 300);

    %%%%%%%%%%
    % Plot for kyy-direction
    kyy_fig = figure;
    h3 = barh(kyy_effect_lower, 'FaceColor', cmap(2, :));
    hold on;
    h4 = barh(kyy_effect_upper, 'FaceColor', cmap(10, :));
    bh3 = get(h3,'BaseLine');
    set(bh3,'BaseValue', kyy_baseline);
    graph_title = ['Sensitivity Analysis (2 layers, 1k simulations)' ' for ' graph_names{g} ' in kyy-direction'];
    title(graph_title, 'FontSize', 13);
    % title('Sensitivity Analysis (2 layers, cs/s, kyy)');

    yticks(1:numel(var_names));
    yticklabels(var_names);
    ax = gca;
    ax.YAxis.FontSize = 14; % for readability purpose
    
    xlim([-7 2.5]);
    xticks([-7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5]);
    xlabel('Mean Permeability');
    ax = gca;
    ax.XAxis.FontSize = 12; % for readability purpose
    ax.XLabel.FontSize = 12;
    legend([h3 h4],'Negative Perturbation','Positive Perturbation')
    % text(h3.XEndPoints, h3.YEndPoints, string(h3.YEndPoints), 'Horiz','right', 'Vert','middle')
    % text(h4.XEndPoints, h4.YEndPoints, string(h4.YEndPoints), 'Horiz','left', 'Vert','middle')
    
    kyy_fig.Position = [100, 100, 800, 400];

    filename = ['abs_' output_name '_kyy_2layer_1k.png'];
    exportgraphics(kyy_fig, filename, 'Resolution', 300);
    % saveas(kyy_fig, filename);

    %%%%%%%%%%
    % Plot for kzz-direction
    kzz_fig = figure;
    h5 = barh(kzz_effect_lower, 'FaceColor', cmap(2, :));
    hold on;
    h6 = barh(kzz_effect_upper, 'FaceColor', cmap(10, :));
    bh5 = get(h5,'BaseLine');
    set(bh5,'BaseValue', kzz_baseline);

    graph_title = ['Sensitivity Analysis (2 layers, 1k simulations)' ' for ' graph_names{g} ' in kzz-direction'];
    title(graph_title, 'FontSize', 13);
    % title('Sensitivity Analysis (2 layers, cs/s, kzz)');

    yticks(1:numel(var_names));
    yticklabels(var_names);
    ax = gca;
    ax.YAxis.FontSize = 14; % for readability purpose

    xlim([-7 2.5]);
    xticks([-7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5]);
    xlabel('Mean Permeability');
    ax = gca;
    ax.XAxis.FontSize = 12; % for readability purpose
    ax.XLabel.FontSize = 12;
    legend([h5 h6],'Negative Perturbation','Positive Perturbation')
    % text(h5.XEndPoints, h5.YEndPoints, string(h5.YEndPoints), 'Horiz','right', 'Vert','middle')
    % text(h6.XEndPoints, h6.YEndPoints, string(h6.YEndPoints), 'Horiz','left', 'Vert','middle')
    
    kzz_fig.Position = [100, 100, 800, 400];
    
    filename = ['abs_' output_name '_kzz_2layer_1k.png'];
    exportgraphics(kzz_fig, filename, 'Resolution', 300);

    close all force;
end