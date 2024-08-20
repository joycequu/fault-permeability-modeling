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

input_data = zeros(11, 5);
input_data(1, :) = [0.2, 0.65, 65, 600, 2000]; % "baseline" 

input_data(2, :) = [0.1, 0.65, 65, 600, 2000]; % "vcl_fw_lower" % sand
input_data(3, :) = [0.3, 0.65, 65, 600, 2000]; % "vcl_fw_upper"

input_data(4, :) = [0.2, 0.4, 65, 600, 2000]; % "vcl_hw_lower" % clay
input_data(5, :) = [0.2, 0.9, 65, 600, 2000]; % "vcl_hw_upper"

input_data(6, :) = [0.2, 0.65, 50, 600, 2000]; % "fdip_lower"
input_data(7, :) = [0.2, 0.65, 80, 600, 2000]; % "fdip_upper"

input_data(8, :) = [0.2, 0.65, 65, 200, 2000]; % "zf_lower"
input_data(9, :) = [0.2, 0.65, 65, 1000, 2000]; % "zf_upper"

input_data(10, :) = [0.2, 0.65, 65, 600, 1000]; % "zmax_lower"
input_data(11, :) = [0.2, 0.65, 65, 600, 3000]; % "zmax_upper"


% n_combo = numel(svcl_fw)*numel(svcl_hw)*numel(sfdip)*numel(szf)*numel(szmax);
% 
% vcl_all_fw = repmat(svcl_fw, n_combo/numel(svcl_fw), 1);
% vcl_all_hw = repmat(svcl_hw, n_combo/numel(svcl_hw), 1);
% fdip_all = repmat(sfdip, n_combo/numel(sfdip), 1);
% zf_all = repmat(szf, n_combo/numel(szf), 1);
% zmax_all = repmat(szmax, n_combo/numel(szmax), 1);
% 
% input_data = [vcl_all_fw vcl_all_hw fdip_all zf_all zmax_all];
% 
% faults_val = nan(100, 6, n_combo);
% faults_grid_isSmear = cell(n_combo, 100);
% faults_grid_units = cell(n_combo, 100);

faults_perm_mean = zeros(11, 3);

for r = 1 : 11
    
    % 2.1 Mandatory input parameters
    % Footwall first and hangingwall next, e.g. {[footwall, FW], [hangingwall, HW]}. 
    % We need to define the layer thickness, clay content, layer dip angle, fault 
    % dip angle, faulting depth, and burial depth. Further details about input parameter 
    % formatting, etc can always be checked from the documentation in the classes 
    % and functions.
    thickness = {100, 100};                % [m]
    vcl       = {input_data(r, 1), input_data(r, 2)};                        % fraction [-]
    dip       = [0, 0];                                                        % [deg.]
    faultDip  = input_data(r, 3);                                                             % [deg.]
    zf        = [input_data(r, 4), input_data(r, 4)];                                                     % [FW, HW], [m]
    zmax      = {input_data(r, 5), input_data(r, 5)};   % {FW, HW}
    dim       = 3;                    % dimensions (2 = 2D, 3 = 3D)
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

        totalPerm = totalPerm + faults{n}.Perm;
        % faults_val(n,:,r) = [myFault.Thick, myFault.Vcl, myFault.Poro, myFault.Perm];
        % faults_grid_units{r, n} = myFault.Grid.units;
        % faults_grid_isSmear{r, n} = myFault.Grid.isSmear;
    end
    faults_perm_mean(r, :) = totalPerm / Nsim;

    telapsed = toc(tstart);
end

% save('threedim_output_two_layers.mat', 'faults_val', 'faults_grid_units', 'faults_grid_isSmear');
save('1layer_faults_perm_mean_output.mat', 'faults_perm_mean');

%% 3. Output Analysis
% 3.1 Visualize stratigraphy and fault (with thickness corresponding to 1st realization)
mySect.plotStrati(faults{1}.Thick, faultDip, unit_plot);  

% 3.2 Visualize intermediate variables
% We define a given parent material (id from 1 to n of materials in stratigraphy), 
% and generate histograms and correlation matrix plots.
layerId = 4;                                            
plotMatPropsHist(faultSections, smears, mySect, layerId, dim) 
% MatProps correlations
[R, P] = plotMatPropsCorr(faultSections, mySect, layerId, dim);
if dim==3
    plotSeg(faults, nSeg)
end

% 3.3 Visualize fault materials
% Visualization for one realization. Choice can be 'randm' (random), 'maxX' 
% (realization with maximum upscaled permeability in across the fault), 'minX', 
% 'maxZ' or 'minZ'.
% General fault materials and perm view
plotId = selectSimId('randm', faults, Nsim);                % simulation index
%plotId = 1;
if U_flex
    faults{plotId}.plotMaterials(faultSections{plotId}{1}, mySect, ...
                                 unit_plot, Us{plotId}) 
else
    faults{plotId}.plotMaterials(faultSections{plotId}{1}, mySect, unit_plot, U) 
end

% 3.4. Visualize upscaled permeability
% Plot upscaled permeability distributions (all simulations)
plotUpscaledPerm(faults, dim)
