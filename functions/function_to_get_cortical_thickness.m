function [thickness, surface, volume, normalised_volume, roiname, global_vals] = function_to_get_cortical_thickness(freesurfer_folder)
% FUNCTION_TO_GET_CORTICAL_THICKNESS Extracts cortical thickness and volume metrics from FreeSurfer stats
%
% Inputs:
%   freesurfer_folder - Path to FreeSurfer subject folder containing stats files
%
% Outputs:
%   thickness - Vector of cortical thickness values (mm) for each ROI
%   surface - Vector of surface area values (mm^2) for each ROI  
%   volume - Vector of volume values (mm^3) for each ROI
%   normalised_volume - Vector of ROI volumes normalized by total cortical volume (%)
%   roiname - Cell array of ROI names
%   global_vals - Structure containing global brain metrics (ICV, cortical vol, WM vol)
%
% Extracts ROI-wise metrics from BB38 atlas stats files and global metrics from
% brainvol.stats. Processes left and right hemispheres separately.

    % Initialize output variables
    thickness = []; surface = []; volume = []; ICV = []; normalized_cortex_gm_volume = [];
    hemispheres = {'lh','rh'};
    
    % Process each hemisphere
    for hem = 1:length(hemispheres)
            clear C
            % Load BB38 atlas stats file
            BB38_stats_file = [freesurfer_folder,'/stats/',hemispheres{hem},'.BB38chimp.stats'];
            if exist(BB38_stats_file) == 2
                fid = fopen(BB38_stats_file);
                % Read ROI metrics (columns: name, surface area, volume, thickness, etc)
                C = textscan(fid,'%s %f %f %f %f %f %f %f %f %f', 'headerlines', 60);
                thickness = [thickness; C{5}]; % Cortical thickness in mm
                surface = [surface; C{3}];     % Surface area in mm^2
                volume = [volume; C{4}];       % Volume in mm^3
                roiname = C{1};
            end
    end
    
    % Load global brain metrics from brainvol.stats
    statsfile = [freesurfer_folder,'/stats/brainvol.stats'];

    if exist(statsfile) == 2
        % Extract total brain volume (ICV)
        [a b] = system(['cat ',statsfile,' | grep "Measure BrainSeg, BrainSegVol, Brain Segmentation Volume"']);
        ss1 =  strsplit(b,'ume,');
        ss2 = strsplit(ss1{2}, ', mm');
        global_vals.ICV = str2double(ss2{1});
        
        % Extract total cortical gray matter volume
        [a b] = system(['cat ',statsfile,' | grep "Measure Cortex, CortexVol, Total cortical gray matter volume"']);
        ss1 =  strsplit(b,'ume,');
        ss2 = strsplit(ss1{2}, ', mm');
        global_vals.cortical_vol = str2double(ss2{1});
        
        % Calculate ROI volumes as percentage of total cortical volume
        normalised_volume = 100 * volume / global_vals.cortical_vol;
        
        % Extract total white matter volume
        [a b] = system(['cat ',statsfile,' | grep "Measure CerebralWhiteMatter, CerebralWhiteMatterVol, Total cerebral white matter volume"']);
        ss1 =  strsplit(b,'ume,');
        ss2 = strsplit(ss1{2}, ', mm');
        global_vals.wm_vol = str2double(ss2{1});
    end