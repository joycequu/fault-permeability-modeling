data = load('fault_perm_3layer_jul7.mat', '-mat');
input = data.all_ee;

% previously: all_ee = zeros(4, 20, 6, 3);
% all_ee = zeros(baseline_selection, 4, 7, 3);

output = zeros(4, 500, 6, 3);

for p1 = 1 : 4
    for p2 = 1 : 500
        for p3 = 1 : 6
            for p4 = 1 : 3
                output(p1, p2, p3, p4) = input(p2, p1, p3, p4);
            end
        end
    end
end

save('fault_perm_3_layer_500.mat', 'output');