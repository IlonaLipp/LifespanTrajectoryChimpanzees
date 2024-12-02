% This script generates surface plots comparing human and chimpanzee brain data
% Key features:
% - Processes multiple contrasts (R1, MTsat, R2s, thickness)
% - Generates statistical maps (z-scores, clusters, etc)
% - Creates visualizations using BrewerMap colormaps
% - Saves results for both hemispheres

clear all
close all
clc

% Add required paths for FreeSurfer MATLAB tools and BrewerMap colormaps
addpath(genpath('freesurfer/6.0.0/ubuntu-xenial-amd64/matlab'))
addpath(genpath('Software/BrewerMap'))

% Define analysis parameters
contrasts = {'R1','MTsat','R2s','thickness'};
hemispheres = {'lh','rh'};

% Process each contrast type
for c = 1:length(contrasts)    
    contrast = contrasts{c};
              
    % Process both hemispheres
    for hem = 1:length(hemispheres)        
        hemisphere = hemispheres{hem};

        % Set cortical depth - use average for thickness, 0.5 for other contrasts
        if strcmp(contrast, 'thickness')     
            cd = 'avg';
        else
            cd = '0.5';
        end
        
        % Set up file paths for results
        foldername = ['/data/pt_02101/results/Cortical_analysis/group_freesurfer/chimps_vs_humans/target_fsaverage/target_fsaverage_',hemisphere,'.',contrast,'_',cd,'_sm03.osgm/humans_vs_chimps'];
        surface_file = ['/data/pt_02101/results/Cortical_analysis/group_freesurfer/chimps_vs_humans/fsaverage/surf/', hemisphere,'.inflated'];
     
        % Copy and rename cluster statistics files
        orig_fname =  [foldername, '/cache.th30.abs.sig.cluster.summary'];
        new_fname =  ['Downloads/cluster_stats_',contrast,'_',hemisphere,'_',cd]
        system(['cp ', orig_fname, ' ', new_fname])
        
        orig_fname =  [foldername, '/cache.th30.abs.sig.cluster.mgh'];
        new_fname =  ['Downloads/Supplementary_Humans_vs_chimpanzees_clusters_',contrast,'_',hemisphere,'.mgh']
        system(['cp ', orig_fname, ' ', new_fname])

        % Generate different types of visualization plots
        vars_to_plot = {'z', 'clusters', 'z_in_clusters', 'chimps_converted_only', 'humans_only'};
        for v = 1:length(vars_to_plot)
            var_to_plot = vars_to_plot{v};
            
            % Load and process data based on plot type
            if strcmp(var_to_plot, 'z')
               [vol_to_plot, M, mr_parms, volsz] = load_mgh([foldername, '/z.mgh']);
               display_range = [-4 4];
            elseif strcmp(var_to_plot, 'clusters')
               [vol_to_plot, M, mr_parms, volsz] = load_mgh([foldername, '/cache.th30.abs.sig.cluster.mgh']);
               display_range = [-4 4];
            elseif strcmp(var_to_plot, 'z_in_clusters')
               [vol_to_plot, M, mr_parms, volsz] = load_mgh([foldername, '/z.mgh']);
               [clu, M, mr_parms, volsz] = load_mgh([foldername, '/cache.th30.abs.sig.cluster.mgh']);
               vol_to_plot(clu == 0) = 0; % Mask non-significant clusters
               display_range = [-4 4];
            elseif strcmp(var_to_plot, 'chimps_converted_only')
               [vol_to_plot, M, mr_parms, volsz] = load_mgh(['Cortical_analysis/group_freesurfer/chimps_vs_humans/target_fsaverage/target_fsaverage_',hemisphere,'.',contrast,'_',cd,'_sm03.osgm/only_chimps/gamma.mgh']);
               display_range = [prctile(vol_to_plot,10) prctile(vol_to_plot,90)];
            elseif strcmp(var_to_plot, 'humans_only')
               [vol_to_plot, M, mr_parms, volsz] = load_mgh(['Cortical_analysis/group_freesurfer/chimps_vs_humans/target_fsaverage/target_fsaverage_',hemisphere,'.',contrast,'_',cd,'_sm03.osgm/only_humans/gamma.mgh']);
               % Use previous display range for humans
            end

            % Generate and save surface plot
            title_name = var_to_plot;
            vol_name = [contrast, ' ', var_to_plot];
            figure_file_to_save_prefix = ['Cortical_analysis/group_freesurfer/Group_plots/Group_figure_brewermap_fsaverage_humans_vs_chimps_',contrast,'_',hemisphere,'_',cd,'_', var_to_plot];
            function_to_create_brewermap(surface_file, vol_to_plot, vol_name, title_name, figure_file_to_save_prefix, display_range);
           
        end            
    end
end
