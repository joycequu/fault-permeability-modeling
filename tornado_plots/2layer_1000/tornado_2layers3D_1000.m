% CURRENT MOST UPDATED CODE FOR TORNADO PLOT

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

% % Values we are evaluating for each variable
% svcl_fw = [0.1; 0.2];
% svcl_hw = [0.4; 0.5];
% % two sets - 0.4 to 0.9 and from 0.1 to 0.9
% sfdip = 50;
% szf = 200;
% szmax = 1000;

% change: want to take log first and then avg instead of avg then log
input = zeros(52, 6);

base_sc_c = [0.2, 0.65, 0.65, 65, 600, 2000];
base_sc_s = [0.2, 0.65, 0.2, 65, 600, 2000];
base_cs_c = [0.65, 0.2, 0.65, 65, 600, 2000];
base_cs_s = [0.65, 0.2, 0.2, 65, 600, 2000];

input(1:13, :) = repmat(base_sc_c, 13, 1);
input(14:26, :) = repmat(base_sc_s, 13, 1);
input(27:39, :) = repmat(base_cs_c, 13, 1);
input(40:52, :) = repmat(base_cs_s, 13, 1);

for i = 1 : 52
    if mod(i, 13) == 8
        input(i, 4) = 50;
        input(i+1, 4) = 80;
    elseif mod(i, 13) == 10
        input(i, 5) = 200;
        input(i+1, 5) = 1000;
    elseif mod(i, 13) == 12
        input(i, 6) = 1000;
        input(i+1, 6) = 3000;
    end
end

% FW: Sand, Clay; HW: Clay
input(2, 1) = 0.1;
input(3, 1) = 0.3;
input(4, 2) = 0.4;
input(5, 2) = 0.9;
input(6, 3) = 0.4;
input(7, 3) = 0.9;

% FW: Sand, Clay; HW: Sand
input(15, 1) = 0.1;
input(16, 1) = 0.3;
input(17, 2) = 0.4;
input(18, 2) = 0.9;
input(19, 3) = 0.1;
input(20, 3) = 0.3;

% FW: Clay, Sand; HW: Clay
input(28, 1) = 0.4;
input(29, 1) = 0.9;
input(30, 2) = 0.1;
input(31, 2) = 0.3;
input(32, 3) = 0.4;
input(33, 3) = 0.9;

% FW: Clay, Sand; HW: Sand
input(41, 1) = 0.4;
input(42, 1) = 0.9;
input(43, 2) = 0.1;
input(44, 2) = 0.3;
input(45, 3) = 0.1;
input(46, 3) = 0.3;

save('2layers_1000_faults_perm_mean_input.mat', 'input');


faults_perm_mean = zeros(52, 3);

for r = 1 : 52
    
    % 2.1 Mandatory input parameters
    % Footwall first and hangingwall next, e.g. {[footwall, FW], [hangingwall, HW]}. 
    % We need to define the layer thickness, clay content, layer dip angle, fault 
    % dip angle, faulting depth, and burial depth. Further details about input parameter 
    % formatting, etc can always be checked from the documentation in the classes 
    % and functions.
    thickness = {[50, 50], 100};                             % fw, hw [m]
    vcl       = {[input(r, 1), input(r, 2)], input(r, 3)};   % fraction [-]
    dip       = [0, 0];                                      % [deg.]
    faultDip  = input(r, 4);                                 % [deg.]
    zf        = [input(r, 5), input(r, 5)];                  % [FW, HW],[m]
    zmax      = {repelem(input(r, 6), numel(vcl{1})),
                 repelem(input(r, 6), numel(vcl{2}))};       % {FW, HW}
    dim       = 3;                                           % dimensions
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
    Nsim              = 1000;       % Number of 3D simulations/realizations

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

        totalPerm = totalPerm + log10(faults{n}.Perm/(milli*darcy));
        % faults_val(n,:,r) = [myFault.Thick, myFault.Vcl, myFault.Poro, myFault.Petarm];
        % faults_grid_units{r, n} = myFault.Grid.units;
        % faults_grid_isSmear{r, n} = myFault.Grid.isSmear;
    end
    faults_perm_mean(r, :) = totalPerm / Nsim;

    telapsed = toc(tstart);
end

% save('threedim_output_two_layers.mat', 'faults_val', 'faults_grid_units', 'faults_grid_isSmear');
save('2layers_1000_faults_perm_mean_output.mat', 'faults_perm_mean');
