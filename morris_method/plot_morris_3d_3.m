data = load('scaled_ee_result.mat', '-mat');
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

% Colors
db = [0, 0, 0.5];       % Dark blue (navy)
mb = [0, 0.5, 1];       % Medium blue
lb = [0.5, 0.75, 1];    % Light blue

dg = [0, 0.5, 0];       % Dark green
mg = [0, 1, 0];         % Medium green
lg = [0.5, 1, 0.5];     % Light green

mo = [1, 0.5, 0];       % Medium organge

cname = {'Sand (FW1), Clay (FW2), Clay (HW)', 'Clay (FW1), Sand (FW2), Clay (HW)', 'Sand (FW1), Clay (FW2), Sand (HW)', 'Clay (FW1), Sand (FW2), Sand (HW)'};
fname = {'sc_c', 'cs_c', 'sc_s', 'cs_s'};
dname = {'kxx', 'kyy', 'kzz'};

label_on = 0; % 0 = off, 1 = on
t = 1; % transparency

abs_rel = 'rel'; % options: 'abs', 'rel'

% distribution: 1 - sc/c, 2 - cs/c, 3 - sc/s, 4 - cs/s
for c = 1 : 4
    % direction: 1 - kxx, 2 - kyy, 3 - kzz
    for d = 1 : 3
        % General Plotting (Everything filled circle)
        x = ee(c, :, d*2-1); % Random x values
        y = ee(c, :, d*2); % Random y values
        labels = {'V_{cl} (FW_1)', 'V_{cl} (FW_2)', 'V_{cl} (HW)', 'f_{dip}', 'z_f', 'z_{max}'}; % Labels for points

        % Plot the scatter plot
        % scatter(x, y, 200, 'filled');

        figure;

        % Vcl (FW1)
        x1 = ee(c, 1, d*2-1);
        y1 = ee(c, 1, d*2);
        scatter(x1, y1, 200, 'x', 'LineWidth', 5, 'MarkerEdgeColor', db, 'MarkerEdgeAlpha', t);
        if label_on
            text(x1 + 0.02, y1 + 0.02, labels{1}, 'FontSize', 12);
        end
        hold on; % Keep the current plot

        % Vcl (FW2)
        x2 = ee(c, 2, d*2-1);
        y2 = ee(c, 2, d*2);
        scatter(x2, y2, 200, 'x', 'LineWidth', 5, 'MarkerEdgeColor', mb, 'MarkerEdgeAlpha', t);
        if label_on
            text(x2 + 0.02, y2 + 0.02, labels{2}, 'FontSize', 12);
        end

        % Vcl (HW)
        x3 = ee(c, 3, d*2-1);
        y3 = ee(c, 3, d*2);
        scatter(x3, y3, 200, 'x', 'LineWidth', 5, 'MarkerEdgeColor', lb, 'MarkerEdgeAlpha', t);
        if label_on
            text(x3 + 0.02, y3 + 0.02, labels{3}, 'FontSize', 12);
        end

        % F dip
        x4 = ee(c, 4, d*2-1);
        y4 = ee(c, 4, d*2);
        % scatter(x4, y4, 150, 'filled', 'MarkerFaceColor', mo, 'MarkerFaceAlpha', 0.5);
        scatter(x4, y4, 200, 'o', 'LineWidth', 2, 'MarkerEdgeColor', mo, 'MarkerFaceAlpha', t);
        if label_on
            text(x4 + 0.02, y4 + 0.02, labels{4}, 'FontSize', 12);
        end

        % Z f
        x5 = ee(c, 5, d*2-1);
        y5 = ee(c, 5, d*2);
        scatter(x5, y5, 200, 'o', 'LineWidth', 2, 'MarkerEdgeColor', dg, 'MarkerFaceAlpha', t);
        if label_on
            text(x5 + 0.02, y5 + 0.02, labels{5}, 'FontSize', 14);
        end

        % Z max
        x6 = ee(c, 6, d*2-1);
        y6 = ee(c, 6, d*2);
        scatter(x6, y6, 200, 'o', 'LineWidth', 2, 'MarkerEdgeColor', lg, 'MarkerFaceAlpha', t);
        if label_on
            text(x6, y6 - 0.02, labels{6}, 'FontSize', 14);
        end

        % Add labels to specific points
        % for i = 1:length(labels)
        %     text(x(i) - 0.2, y(i) - 0.2, labels{i}, 'FontSize', 16);
        % end

        % for i = 1:length(labels)
        %     text(x(i), y(i), labels{i}, 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        % end

        axis tight;

        hold off; % Release the plot hold
        xlabel('Mean', 'FontSize', 16);
        ylabel('Standard Deviation','FontSize', 16);
        graph_title = [cname{c} ' - ' dname{d} ' direction'];
        title(graph_title,'FontSize', 16);
        
        if strcmp(abs_rel, 'rel')
            xlim([0, max(x) + (max(x) - min(x)) * 0.1]);
            ylim([0, max(y) + (max(y) - min(y)) * 0.1]);
        else
            xlim([min_mean, max_mean + (max_mean-min_mean)*0.05]);
            ylim([min_std, max_std + (max_std - min_std)*0.05]);
        end

        set(gca, 'FontSize', 12);

        grid on;

        % filename = ['rel_' output_name '_kzz_2layer_1k.png'];
        filename = ['scaled_morris_' abs_rel '_' fname{c} '_' dname{d} '_2layer_500.png'];
        exportgraphics(gcf, filename, 'Resolution', 300);

        % gname(labels);
        clf;
    end
end

close all force;
% Save the plot as PNG
% saveas(gcf, 'scc_kxx.png');