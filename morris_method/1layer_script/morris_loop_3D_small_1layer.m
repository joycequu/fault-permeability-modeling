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
mrstModule add mrst-gui coarsegrid upscaling incomp mpfa mimetic
mrstVerbose on     % set to on for more insight in the command window

%% 2. Define Model and Upscale Permeability

all_ee = zeros(1000, 5, 3);
% 5 for the parameters
% 20 for 20 runs
% 3 for kxx, kyy, kzz

% 1 - vcl_fw
% 2 - vcl_hw
% 3 - sfdip
% 4 - szf
% 5 - szmax

% change we will implement on each variable
delta = [0, 0, 0, 0, 0;
         0.01, 0, 0, 0, 0;
         0, 0.05, 0, 0, 0;
         0, 0, 0.5, 0, 0;
         0, 0, 0, 2, 0;
         0, 0, 0, 0, 10];

% faults{i}.Grid.units
units = cell(1000, 1);
% perm will be an array
perm = zeros(1000, 3);

% do it 5 times first
for l = 1 : 1000 % l = EEk (number of repeats)
    
    % initialize pbase value to set later on
    pbase = [0, 0, 0];
    
    % generate new baseline value only when it's a new repeat
    fw_val = 0.1 + (0.3 - 0.1)*rand;
    hw_val = 0.5 + (0.9 - 0.5)*rand;
    sfdip_val = 50 + (80 - 50)*rand;
    szf_val = 200 + (1000 - 200)*rand;
    szmax_val = 1000 + (3000 - 1000)*rand;

    for d = 1 : 6

        % faults{i}.Grid.units
        units = cell(1000, 1);
        % perm will be an array
        perm = zeros(1000, 3);
        
        faultPerm = [0, 0, 0];

        % 2.1 Mandatory input parameters
        % Footwall first and hangingwall next, e.g. {[footwall, FW], [hangingwall, HW]}. 
        % We need to define the layer thickness, clay content, layer dip angle, fault 
        % dip angle, faulting depth, and burial depth. Further details about input parameter 
        % formatting, etc can always be checked from the documentation in the classes 
        % and functions.
        thickness = {100, 100};                % [m]
        vcl       = {fw_val + delta(d, 1), hw_val + delta(d, 2)};         % fraction [-]
        dip       = [0, 0];                                               % [deg.]
        faultDip  = sfdip_val + delta(d, 3);                              % [deg.]
        zf        = [szf_val + delta(d, 4), szf_val + delta(d, 4)];       % [FW, HW], [m]
        zmax      = {szmax_val + delta(d, 5), szmax_val + delta(d, 5)};   % {FW, HW}
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
        Nsim              = 100;       % Number of 3D simulations/realizations

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
        faults = cell(Nsim, 1);
        Us = cell(Nsim, 1);
        nSeg_fcn = nSeg.fcn;
        U_flex = U.flexible;
        tstart = tic;
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
            for k=1:numel(myFault.SegLen)
                % Get material property (intermediate variable) samples, and fix
                % fault thickness of current realization (3D only).
                myFaultSection = myFaultSection.getMaterialProperties(mySect, 'corrCoef', rho);
                myFaultSection.MatProps.thick = myFault.Thick;
                if isempty(G)
                    G = makeFaultGrid(myFault.Thick, myFault.Disp, ...
                              myFault.Length, myFault.SegLen, U);
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
            end

            % Compute 3D upscaled permeability distribution
            if U_flex
                myFault = myFault.upscaleProps(G, Us{n});
            else
                myFault = myFault.upscaleProps(G, U);
            end
    
            % Save results
            faults{n} = myFault;
            if mod(n, 50) == 0
                disp(['Simulation ' num2str(n) ' / ' num2str(Nsim) ' completed.'])
            end

            faultPerm = faultPerm + faults{n}.Perm;

        end
        % end of the Nsim



        % CHANGE THIS
        faultPerm = faultPerm / Nsim;
        if d == 1
            pbase = faultPerm;
        else
            % storing the 
            faultPerm = (faultPerm - pbase)/delta(d, d-1);
            all_ee(l, d-1, :) = faultPerm;
        end

        telapsed = toc(tstart);
    end
    % end of loop for different parameters for each run
end

% save('threedim_output_two_layers.mat', 'faults_val', 'faults_grid_units', 'faults_grid_isSmear');
save('all_ee_twenty.mat', 'all_ee');
