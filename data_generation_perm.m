%% Example 0: Single stratigraphic case + analysis (3D)
% This is a complete introductory example. It shows how to load the appropriate 
% MRST modules, define the inputs according to a given faulted stratigraphy, and 
% generate the output permeability distributions. A comprehensive analysis of 
% the results is also shown. The algorithm is run in 3D mode.
% 
% We first make sure that the workspace is clean:
clear
close all force
%rng('default')          % check repeatability

%% 1. Load Required MRST Modules
% First, navigate to the mrst folder and run |startup.m|. We can then load the 
% appropriate modules for generating MRST grids and upscale the permeability:
mrstModule add mrst-gui coarsegrid upscalingincomp mpfa mimetic
mrstVerbose on     % set to on for more insight in the command window

%% 2. Define Model and Upscale Permeability

baseline_selection = 500;

% all_ee = zeros(baseline_selection, 4, 7, 3);
% 4 for the sand/clay combos
% 6 for the parameters
% 20 for 20 runs
% 3 for kxx, kyy, kzz

% 1 - vcl_fw1
% 2 - vcl_fw2
% 3 - vcl_hw
% 4 - fdip
% 5 - zf
% 6 - zmax

% repeat four times for each different combination of 3layer distribution

% need to modify this to work for 3 layers

% change we will implement on each variable

seq_num = 1;

delta = [0, 0, 0, 0, 0, 0;
         0.02, 0, 0, 0, 0, 0;
         0, 0.04, 0, 0, 0, 0;
         0, 0, 0.04, 0, 0, 0;
         0, 0, 0, 3, 0, 0;
         0, 0, 0, 0, 70, 0;
         0, 0, 0, 0, 0, 180];

for combo = 1 : 4
    if combo == 1 % sand (fw1), clay (fw2), clay (hw)
        vcl_val = [0.1, 0.5, 0.5;
                   0.2, 0.4, 0.4];
    elseif combo == 2 % clay (fw1), sand (fw2), clay (hw)
        delta(2, 1) = 0.04;
        delta(3, 2) = 0.02;
        delta(4, 3) = 0.04;
        vcl_val = [0.5, 0.1, 0.5;
                   0.4, 0.2, 0.4];
    elseif combo == 3 % sand (fw1), clay (fw2), sand (hw)
        delta(2, 1) = 0.02;
        delta(3, 2) = 0.04;
        delta(4, 3) = 0.02;
        vcl_val = [0.1, 0.5, 0.1;
                   0.2, 0.4, 0.2];
    else % clay (fw1), sand (fw2), sand (hw)
        delta(2, 1) = 0.04;
        delta(3, 2) = 0.02;
        delta(4, 3) = 0.02;
        vcl_val = [0.5, 0.1, 0.1;
                   0.4, 0.2, 0.2];
    end

    for l = 1 : baseline_selection % l = EEk (number of repeats)
    
        % initialize pbase value to set later on
        pbase = [0, 0, 0];
    
        % generate new baseline value only when it's a new repeat
        fw_val1 = vcl_val(1, 1) + vcl_val(2, 1)*rand;
        fw_val2 = vcl_val(1, 2) + vcl_val(2, 2)*rand;
        hw_val = vcl_val(1, 3) + vcl_val(2, 3)*rand;
        fdip_val = 50 + (80 - 50)*rand;
        zf_val = 200 + (900 - 200)*rand;
        zmax_val = 1000 + (2800 - 1000)*rand;

        for d = 1 : 7

            % 2.1 Mandatory input parameters
            thickness = {[50, 50], 100};                  % [m]
            vcl       = {[fw_val1 + delta(d, 1), fw_val2 + delta(d, 2)], hw_val + delta(d, 3)};         % fraction [-]
            dip       = [0, 0];                                               % [deg.]
            faultDip  = fdip_val + delta(d, 4);                              % [deg.]
            zf        = [zf_val + delta(d, 5), zf_val + delta(d, 5)];       % [FW, HW], [m]
            zmax      = {repelem(zmax_val + delta(d, 6), numel(vcl{1})), repelem(zmax_val + delta(d, 6), numel(vcl{2}))};   % {FW, HW}
            dim       = 3;                                                    % dimensions (2 = 2D, 3 = 3D)
            unit_plot = 'm';

            % 2.2 Optional input parameters
            % In this case, we indicate a maximum fault material permeability of and a correlation 
            % coefficient for dependent variables:
            maxPerm = 1000;                 % [mD]
            rho     = 0.6;                  % Corr. coeff. for multivariate distributions
            totalPerm = 0;

            % 2.3 Flow upscaling options and number of simulations
            U.method          = 'tpfa';     % 'tpfa'for 3D
            U.coarseDims      = [1 1 1];    % desired n cells [x, y, z] in coarse grid
            U.flexible        = true;       % default true, much faster but U.coarseDims
                                            % will be modified in some realizations
                                            % unless U.coarseDims = [1 1 1] (do not set
                                            % to false in that case).
            Nsim              = 200;       % Number of 3D simulations/realizations
            
            % 2.4 Define Stratigraphy and FaultedSection objects
            % Organize the input parameters in HW and FW, and use that info to create a 
            % FaultedSection object which contains all required information.
            % FW and HW
            footwall = Stratigraphy(thickness{1}, vcl{1}, 'Dip', dip(1), ...
                            'DepthFaulting', zf(1), 'DepthBurial', zmax{1});
            hangingwall = Stratigraphy(thickness{2}, vcl{2}, 'Dip', dip(2), 'IsHW', 1, ...
                            'NumLayersFW', footwall.NumLayers, ...
                            'DepthFaulting', zf(2), 'DepthBurial', zmax{2});

            % Instantiate FaultedSection object (Strati in Faulted Section)
            mySect = FaultedSection(footwall, hangingwall, faultDip, 'maxPerm', maxPerm);

            % 2.5 Get material distributions
            % We use the inputs to constrain the ranges and distributions for each of the 
            % intermediate variables.
            % Get material property distributions
            mySect = mySect.getMatPropDistr();
            % Get along-strike segmentation
            nSeg = getNSeg(mySect.Vcl, mySect.IsClayVcl, mySect.DepthFaulting);

            % 2.6 Generate intermediate variable samples, calculate smear dimensions 
            %     and upscale permeability.
            % We create two container variables (faults and smears) where we'll save all 
            % data for each realization. For each realization, the code defines a Fault object, 
            % generates intermediate variable samples, calculates the smear dimensions, and, 
            % within upscaleSmearPerm, generates a fault material distribution consistent 
            % with the inputs and upscales the permeability.
            % Generate fault object with properties for each realization
            assert(dim==3);
            faultSections = cell(Nsim, 1);
            smears = cell(Nsim, 1);
            Us = cell(Nsim, 1);
            nSeg_fcn = nSeg.fcn;
            U_flex = U.flexible;
            tstart = tic;

            % Need to save output for each variation of input
            perm_data = cell(Nsim, 4);
            units_data = cell(Nsim, 1);
            % upscale permeability saved for sorting purposes
            upscale_perm = zeros(Nsim, 3);

            parfor n=1:Nsim    % parfor allowed if you have the parallel computing toolbox
            %for n=1
            % Instantiate fault section and get segmentation for this realization
                myFaultSection = Fault2D(mySect, faultDip);
                myFault = Fault3D(myFaultSection, mySect);
                if U_flex
                    [myFault, Us{n}] = myFault.getSegmentationLength(U, nSeg_fcn);
                else
                    myFault = myFault.getSegmentationLength(U, nSeg_fcn);
                end
                G = [];

                % temp_perm to resolve dependency issue
                temp_perm = cell(1, 4);

                for k=1:numel(myFault.SegLen)
                    % Get material property (intermediate variable) samples, and fix
                    % fault thickness of current realization (3D only).
                    myFaultSection = myFaultSection.getMaterialProperties(mySect, 'corrCoef', rho);
                    myFaultSection.MatProps.thick = myFault.Thick;
                    if isempty(G)
                        G = makeFaultGrid(myFault.Thick, myFault.Disp, ...
                                myFault.Length, myFault.SegLen, U);
                        [nx, nz] = deal(G.cartDims(1), G.cartDims(3));
                        % reformatted
                        for w = 1:4
                            temp_perm{w} = zeros(nx, nz, 16);
                        end
                        % [perm_data{n,:}] = deal(zeros(nx, nz, 16)); % preallocate each 3D matrix for a given realization n
                    end
        
                    % Generate smear object with T, Tap, L, Lmax
                    smear = Smear(mySect, myFaultSection, G, 1);
        
                    % Place fault materials and assign cell-based properties in 2D section
                    myFaultSection = myFaultSection.placeMaterials(mySect, smear, G);
        
                    % Extrude 2D section to fill current segment
                    myFault = myFault.assignExtrudedVals(G, myFaultSection, k);
        
                    % Save results
                    faultSections{n}{k} = myFaultSection;
                    smears{n}{k} = smear;
                    
                    temp_perm{1}(:, :, k) = reshape(myFaultSection.Grid.perm(:,1), nx, nz); % kxx
                    temp_perm{2}(:, :, k) = reshape(myFaultSection.Grid.permy, nx, nz); % kyy
                    temp_perm{3}(:, :, k) = reshape(myFaultSection.Grid.perm(:,3), nx, nz); % kzz
                    temp_perm{4}(:, :, k) = reshape(myFaultSection.Grid.perm(:,2), nx, nz); % kxz
 
                end

                % Compute 3D upscaled permeability distribution
                if U_flex
                    myFault = myFault.upscaleProps(G, Us{n});
                else
                    myFault = myFault.upscaleProps(G, U);
                end
    
                if mod(n, 50) == 0
                    disp(['Simulation ' num2str(n) ' / ' num2str(Nsim) ' completed.'])
                end

                % add up the fault permeability for morris method analysis
                % faultPerm = faultPerm + log10(myFault.Perm/(milli*darcy));

                % save units & perm for data generation purposes
                % units{n} = myFault.Grid.units;
                % perms{n} = myFault.Grid.perms;
                % perm(n, :) = myFault.Perm;

                upscale_perm(n, :) = log10(myFault.Perm/(milli*darcy));
                perm_data(n, :) = temp_perm;
                units_data{n} = myFault.Grid.units;

            end
            % end of the Nsim

            upscale_perm_indexed = [(1:Nsim)', zeros(Nsim, 1), zeros(Nsim, 1), zeros(Nsim, 1), zeros(Nsim, 1)];
            upscale_perm_index(:, 2:4) = upscale_perm;
            
            % 9 - kxx - P10, Median, P90; kyy - ...; kzz - ...
            perm_output = cell(9, 4);
            units_output = cell(9, 1);
            upscale_perm_output = zeros(9, 3);

            for dim = 2 : 4
                sorted_upscale_perm = sortrows(upscale_perm_indexed, dim);
                p10_index = sorted_upscale_perm(Nsim * 10, 1);
                p50_index = sorted_upscale_perm(Nsim * 50, 1);
                p90_index = sorted_upscale_perm(Nsim * 90, 1);
                % save to permeability output
                perm_output((dim-2)*3 + 1) = perm_data(p10_index, :);
                perm_output((dim-2)*3 + 2) = perm_data(p50_index, :);
                perm_output((dim-2)*3 + 3) = perm_data(p90_index, :);
                % save to units output
                units_output((dim-2)*3 + 1) = units_data{p10_index};
                units_output((dim-2)*3 + 2) = units_data{p50_index};
                units_output((dim-2)*3 + 3) = units_data{p90_index};
                % save to upscale permeability output
                upscale_perm_output((dim-2)*3 + 1, :) = upscale_perm(p10_index, :);
                upscale_perm_output((dim-2)*3 + 2, :) = upscale_perm(p50_index, :);
                upscale_perm_output((dim-2)*3 + 3, :) = upscale_perm(p90_index, :);
            end

            inputs = cell(6, 1);
            inputs{1} = thickness;
            inputs{2} = vcl;
            inputs{3} = dip;
            inputs{4} = faultDip;
            inputs{5} = zf;
            inputs{6} = zmax;

            telapsed = toc(tstart);

            filename = strcat("data_3layer_seq_full_", num2str(seq_num), ".mat");
            seq_num = seq_num + 1;

            save(filename, 'inputs', 'perm_output', 'units_output', 'upscale_perm_output'); % removed 'units'

            clear nx nz myFault perm_output units_output upscale_perm_output upscale_perm upscale_perm_index perm_data units_data inputs
        end
        % end of loop for different parameters for each run
    end
end

% save('threedim_output_two_layers.mat', 'faults_val', 'faults_grid_units', 'faults_grid_isSmear');
% save('scaled_fault_perm_3_layer.mat', 'all_ee');