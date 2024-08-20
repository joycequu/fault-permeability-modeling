function plotUpscaledPerm(faults, dim, plotOpt, fignum)
%
%
%

%% Preliminaries

% Utilities
latx = {'Interpreter', 'latex'};
sz = [14, 12];

% Fault MatProps
if nargin < 2 || dim == 2
    kAlongStrike = cell2mat(cellfun(@(x) x.Grid.permy, faults, ...
                                    'UniformOutput', false)')./(milli*darcy);
end
perms = cell2mat(cellfun(@(x) x.Perm, faults, ...
                         'UniformOutput', false)) ./ (milli*darcy);
s = size(faults{1}.Perm, 1);                     
if s > 1
    perm_cell = zeros(numel(faults), dim, s);
    for n=1:s
        perm_cell(:,:,n) = cell2mat(cellfun(@(x) x.Perm(n,:), faults, ...
                                 'UniformOutput', false)) ./ (milli*darcy);
    end
    K_cell = log10(perm_cell);
end
if nargin < 2 || dim == 2
    thick = cell2mat(cellfun(@(x) x.MatProps.thick, faults, ...
                             'UniformOutput', false));
elseif dim == 3
    thick = cell2mat(cellfun(@(x) x.Thick, faults, 'UniformOutput', false));
end
vcl = cell2mat(cellfun(@(x) x.Vcl, faults, 'UniformOutput', false));
if any(any(perms < 0))
    id = unique([find(perms(:, 1)<0), find(perms(:, 2)<0), find(perms(:, 3)<0)]);
    warning(['Negative upscaled perms found in ' num2str(numel(id))...
             ' simulations (ignored).'])
    kAlongStrike(:, id) = [];
    perms(id, :) = [];
    thick(id, :) = [];
    vcl(id, :) = [];
end
if nargin < 2 || dim == 2
    logkStrikeBounds = log10(getAveragingPerm(kAlongStrike, {'ha', 'ah'}));
end
    
% Hist params
K = log10(perms);
nbins = 25;
logMinP = -6; %min(min(K));
logMaxP = 2; %max(max(K));
edges = linspace(fix(logMinP)-1, fix(logMaxP)+1, nbins);

% Plot
if nargin > 2 && strcmp(plotOpt, 'histOnly')
    % Histograms
    if nargin < 4
        fh = figure(randi(10000, 1, 1));
        %tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
        nexttile(1)
    else
        fh = figure(fignum);
        hold on
    end
    labls = ["$\log_{10}(k_{xx}$ [mD])", ...
        "$\log_{10}(k_{yy}$ [mD])", ...
        "$\log_{10}(k_{zz}$ [mD])"];
    histogram(K(:, 1), edges, 'Normalization', 'probability', ...
        'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
    xlabel(labls(1), latx{:}, 'fontSize', sz(2))
    ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    xlim([fix(logMinP)-1 fix(logMaxP)+1])
    ylim([0 1]); yticks(0:.2:1)
    grid on
    %xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
    
    nexttile(2)
    rr = [255, 125, 125]/255;
    if nargin < 2 || dim == 2
    plot(repelem(logkStrikeBounds(1), 2), [0 1], '-', ...
         'color', 'k', 'lineWidth', 1);
    hold on
    plot(logkStrikeBounds(1), 0.5, 'dk', 'markerFacecolor', rr, 'markerSize', 4)
    plot(repelem(logkStrikeBounds(2), 2), [0 1], '-', ...
         'color', 'k', 'lineWidth', 1);
    plot(logkStrikeBounds(2), 0.5, 'dk', 'markerFacecolor', rr, 'markerSize', 4)
    end
    histogram(K(:, 2), edges, 'Normalization', 'probability', ...
        'FaceColor', rr, 'FaceAlpha', 1)
    xlabel(labls(2), latx{:}, 'fontSize', sz(2))
    %ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    xlim([fix(logMinP)-1 fix(logMaxP)+1])
    ylim([0 1]); yticks(0:.2:1)
    grid on
    %xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
    hold off
    
    nexttile(3)
    bb = [125, 125, 255]/255;
    histogram(K(:, 3), edges, 'Normalization', 'probability', ...
        'FaceColor', bb, 'FaceAlpha', 1)
    xlabel(labls(3), latx{:}, 'fontSize', sz(2))
    %ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    xlim([fix(logMinP)-1 fix(logMaxP)+1])
    ylim([0 1]); yticks(0:.2:1)
    grid on
    %xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
    %set(fh, 'position', [200, 200, 150, 350]);
    set(fh, 'position', [200, 200, 600, 200]);
    
    if s > 1 && s <= 10
        fh = figure(randi(10000, 1, 1));
        tiledlayout(s, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
        for n=1:s
            nexttile(1+3*(n-1))
            labls = ["$\log_{10}(k_{xx}$ [mD])", ...
                "$\log_{10}(k_{yy}$ [mD])", ...
                "$\log_{10}(k_{zz}$ [mD])"];
            histogram(K_cell(:, 1, n), edges, 'Normalization', 'probability', ...
                'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
            if n==1
                xlabel(labls(1), latx{:}, 'fontSize', sz(2))
                ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            end
            xlim([fix(logMinP)-1 fix(logMaxP)+1])
            ylim([0 1]); yticks(0:.2:1)
            grid on
            %xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
            
            nexttile(2+3*(n-1))
            histogram(K_cell(:, 2, n), edges, 'Normalization', 'probability', ...
                'FaceColor', rr, 'FaceAlpha', 1)
            if n == 1
                xlabel(labls(2), latx{:}, 'fontSize', sz(2))
            end
            %ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            xlim([fix(logMinP)-1 fix(logMaxP)+1])
            ylim([0 1]); yticks(0:.2:1)
            grid on
            %xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
            hold off
            
            nexttile(3+3*(n-1))
            histogram(K_cell(:, 3, n), edges, 'Normalization', 'probability', ...
                'FaceColor', bb, 'FaceAlpha', 1)
            if n==1
                xlabel(labls(3), latx{:}, 'fontSize', sz(2))
            end
            %ylabel('P [-]', latx{:}, 'fontSize', sz(2))
            xlim([fix(logMinP)-1 fix(logMaxP)+1])
            ylim([0 1]); yticks(0:.2:1)
            grid on
            %xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
            %set(fh, 'position', [200, 200, 150, 350]);
        end
        set(fh, 'position', [100, 100, 600, 1200]);
    end
    
else
    if nargin > 2 && strcmp(plotOpt, 'all')
        % fault thickness  vs perm
        [Rx, Px] = corr(K(:,1), thick, 'Type', 'Spearman');         % corrcoeff and pval matrices
        [Ry, Py] = corr(K(:,2), thick, 'Type', 'Spearman');
        [Rz, Pz] = corr(K(:,3), thick, 'Type', 'Spearman');
        a = 0.05;                                    % significance level
        pvals = [Px Py Pz];
        r = [Rx Ry Rz];
        if pvals(1) < a, colr{1} = 'k'; fw{1} = 'bold';
        else, colr{1} = 'k'; fw{1} = 'normal'; end
        if pvals(2) < a, colr{2} = 'r'; fw{2} = 'bold';
        else, colr{2} = 'r'; fw{2} = 'normal'; end
        if pvals(3) < a, colr{3} = 'b'; fw{3} = 'bold';
        else, colr{3} = 'b'; fw{3} = 'normal'; end
        fh = figure(randi(10000, 1, 1));
        tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        nexttile(1)
        hold on
        scatter(thick, K(:,1), 8, 'ok','MarkerEdgeAlpha', 0.6)
        scatter(thick, K(:,2), 8, 'sr','MarkerEdgeAlpha', 0.3)
        scatter(thick, K(:,3), 8, 'db','MarkerEdgeAlpha', 0.15)
%         if any(pvals ~= 0) && any(log10(pvals) > -10)
%             title({['$r_s = $ ' num2str(round(r(1), 3)) ', ' ...
%                    num2str(round(r(2), 3)) ', ' num2str(round(r(3), 3))]; ...
%                    ['$\log_{10} (p) = $' num2str(round(log10(pvals(1)))) ', ' ...
%                    num2str(round(log10(pvals(2)))) ', ' ...
%                    num2str(round(log10(pvals(3))))]}, ...
%                    latx{:}, 'fontSize', 10);
%         else
            title({['$r_s = $ ' num2str(round(r(1), 3)) ', ' ...
                   num2str(round(r(2), 3)) ', ' num2str(round(r(3), 3))]}, ...
                   latx{:}, 'fontSize', 10);
        %end
        xlabel('f$_\mathrm{T}$ [m]', latx{:}, 'fontSize', sz(2))
        ylabel('$\log_{10}(k$ [mD])', latx{:}, 'fontSize', sz(2))
        ylim([fix(logMinP)-1 fix(logMaxP)+1])
        xlim([0 fix(max(thick))+1])
        h = legend({'$k_{xx}$', '$k_{yy}$', '$k_{zz}$'}, latx{:}, 'fontSize', 10);
        set(h.BoxFace, 'ColorType','truecoloralpha', ...
            'ColorData', uint8(255*[1;1;1;.9])); 
        grid on
        hold off
        
        % N smear domains vs perm
        [Rx, Px] = corr(K(:,1), vcl, 'Type', 'Spearman');         % corrcoeff and pval matrices
        [Ry, Py] = corr(K(:,2), vcl, 'Type', 'Spearman');
        [Rz, Pz] = corr(K(:,3), vcl, 'Type', 'Spearman');
        a = 0.05;                                    % significance level
        pvals = [Px Py Pz];
        r = [Rx Ry Rz];
        if pvals(1) < a, colr{1} = 'm'; fw{1} = 'bold';
        else, colr{1} = 'k'; fw{1} = 'normal'; end
        if pvals(2) < a, colr{2} = 'm'; 
        else, colr{2} = 'r'; fw{2} = 'normal'; end
        if pvals(3) < a, colr{3} = 'm';
        else, colr{3} = 'b'; fw{3} = 'normal'; end
        nexttile(2)
        hold on
        scatter(vcl, K(:,1), 8, 'ok','MarkerEdgeAlpha', 0.6)
        scatter(vcl, K(:,2), 8, 'sr','MarkerEdgeAlpha', 0.3)
        scatter(vcl, K(:,3), 8, 'db','MarkerEdgeAlpha', 0.15)
%         if any(pvals ~= 0) && any(log10(pvals) > -10)
%             title({['$r_s = $ ' num2str(round(r(1), 3)) ', ' ...
%                    num2str(round(r(2), 3)) ', ' num2str(round(r(3), 3))]; ...
%                    ['$\log_{10} (p) = $' num2str(round(log10(pvals(1)))) ', ' ...
%                    num2str(round(log10(pvals(2)))) ', ' ...
%                    num2str(round(log10(pvals(3))))]}, ...
%                    latx{:}, 'fontSize', 10);
%         else
            title({['$r_s = $ ' num2str(round(r(1), 3)) ', ' ...
                   num2str(round(r(2), 3)) ', ' num2str(round(r(3), 3))]}, ...
                   latx{:}, 'fontSize', 10);
%        end
        xlabel('f$_{V_\mathrm{cl}}$', latx{:}, 'fontSize', sz(2))
        ylabel('$\log_{10}(k$ [mD])', latx{:}, 'fontSize', sz(2))
        ylim([fix(logMinP)-1 fix(logMaxP)+1])
        xlim([min(vcl) max(vcl)])
        %xticks([0 0.2 0.4 0.6 0.8 1])
        grid on
        hold off
        set(fh, 'position', [200, 200, 200, 400]);
    end
    
    % Histograms
    fh = figure(randi(10000, 1, 1));
    tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    labls = ["$\log_{10}(k_{xx}$ [mD])", ...
             "$\log_{10}(k_{yy}$ [mD])", ...
             "$\log_{10}(k_{zz}$ [mD])"];
    nexttile(1)
    histogram(K(:, 1), edges, 'Normalization', 'probability', ...
        'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
    xlabel(labls(1), latx{:}, 'fontSize', sz(2))
    ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    xlim([fix(logMinP)-1 fix(logMaxP)+1])
    ylim([0 1]); yticks(0:.2:1)
    grid on
    %xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
    
    nexttile(2)
    rr = [255, 125, 125]/255;
    if nargin < 2 || dim == 2
        plot(repelem(logkStrikeBounds(1), 2), [0 1], '-', ...
             'color', 'k', 'lineWidth', 1);
        hold on
        plot(logkStrikeBounds(1), 0.3, 'dk', 'markerFacecolor', rr, 'markerSize', 4)
        plot(repelem(logkStrikeBounds(2), 2), [0 1], '-', ...
             'color', 'k', 'lineWidth', 1);
        plot(logkStrikeBounds(2), 0.3, 'dk', 'markerFacecolor', rr, 'markerSize', 4)
    end
    histogram(K(:, 2), edges, 'Normalization', 'probability', ...
        'FaceColor', rr, 'FaceAlpha', 1)
    xlabel('$\hat{k}_{yy}$ [mD]', latx{:}, 'fontSize', sz(2))
    ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    xlim([fix(logMinP)-1 fix(logMaxP)+1])
    ylim([0 1]); yticks(0:.2:1)
    grid on
    %xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
    hold off
    
    nexttile(3)
    bb = [125, 125, 255]/255;
    histogram(K(:, 3), edges, 'Normalization', 'probability', ...
        'FaceColor', bb, 'FaceAlpha', 1)
    xlabel(labls(3), latx{:}, 'fontSize', sz(2))
    ylabel('P [-]', latx{:}, 'fontSize', sz(2))
    xlim([fix(logMinP)-1 fix(logMaxP)+1])
    ylim([0 1]); yticks(0:.2:1)
    grid on
    %xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
    
    % Scatters
    %tidss = [2 3 4 6 7 8];
    tidss = [4 5 6];
    x = [1 2 3];
    y = [2 3 1];
    %colrs = [128, 0, 0; 0, 119, 128; 128, 0, 0; 128, 119, 128; 0, 119, 128; ...
    %    128, 119, 128]./255;
    [R, P] = corrcoef(K);           % corrcoeff and pval matrices
    a = 0.05;                       % significance level
    pvals = P'; pvals = pvals([4,8,3]);
    r = R'; r = r([4, 8, 3]);
    markr = ['+', 'x', '^'];
    for n=1:numel(tidss)
        nexttile(tidss(n))
        colormap(turbo);
%        scatter(K(:,x(n)), K(:,y(n)), 4, 'k.', 'MarkerEdgeAlpha', 0.2)
        histogram2(K(:,x(n)), K(:,y(n)), edges, edges, ...
            'Normalization', 'Probability', 'DisplayStyle','tile', ...
            'ShowEmptyBins','off');
        c = colorbar;
        if n == 1
            c.Label.Interpreter = 'latex';
            c.Label.String = 'P [-]';
            c.Label.FontSize = 9;
        end
        if pvals(n) < a
            colr = 'm';
            fw = 'bold';
        else
            colr = 'k';
            fw = 'normal';
        end
        text(edges(1) + 0.05*(edges(end) - edges(1)), ...
            edges(end) - 0.05*(edges(end) - edges(1)), ...
            ['\rho = ' num2str(round(r(n), 3))],  'color', colr, 'fontSize', 10, 'fontWeight', fw);
        xlim([fix(logMinP)-1 fix(logMaxP)+1])
        ylim([fix(logMinP)-1 fix(logMaxP)+1])
        grid on
        xlabel(labls{x(n)}, latx{:}, 'fontSize', sz(2))
        ylabel(labls{y(n)}, latx{:}, 'fontSize', sz(2))
    end
    %set(fh, 'position', [200, 200, 500, 300]);
    set(fh, 'position', [200, 200, 1400, 700]);
end
end