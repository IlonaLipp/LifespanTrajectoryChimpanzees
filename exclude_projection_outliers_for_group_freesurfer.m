% This script processes MRI data for a group of subjects using FreeSurfer.
% It identifies and corrects outliers in cortical surface projections for 
% different contrasts and depths. The script generates corrected .mgh files 
% and visualizes both uncorrected and corrected data, saving the plots for 
% each subject and contrast.

% Clear workspace and figures
close all
clear all
clc

% Add required FreeSurfer MATLAB tools and BrewerMap colormaps
addpath(genpath('freesurfer/6.0.0/ubuntu-xenial-amd64/matlab'))
addpath(genpath('Software/BrewerMap'))

% Set up analysis parameters
analysis_folder = 'Cortical_analysis/group_freesurfer';
contrasts = {'R1','MTsat','R2s'};
contrast_names = {'R1 (s^{-1}) 7T','MT_{sat} (p.u.) 7T', 'R2* (s^{-1}) 7T''};

% Define display ranges for each contrast
display_ranges = {[1.4 2],[3.5 4.5],[15 20]};
hemispheres = {'lh','rh'};

% Create directory for individual subject plots
plot_folder = [analysis_folder,'/individual_plots'];
system(['mkdir ',plot_folder]);

% Get list of subject directories
subj_dirs = dir([analysis_folder,'/Subj*']);

% Process each subject
for s = 1:length(subj_dirs)
    subj_dir = [subj_dirs(s).folder,'/',subj_dirs(s).name];  

    % Process each contrast type (R1, MTsat, R2s)
    for c = 1:length(contrasts)
      contrast = contrasts{c};
      display_range = display_ranges{c};
      
      % Process different cortical depths (0=average, 0.2-0.8=specific depths)
      for depth = [0 0.2 0.35 0.5 0.65 0.8]
          dstr = num2str(depth);
          
          % Process left and right hemispheres
          for h = 1:2
            hemisphere = hemispheres{h};
            
            % Define file paths for original, corrected and backup files
            mgh_file = [subj_dir,'/Profiles/equi_',hemisphere,'_whole_brain_',contrast,'_0p3_',dstr,'.mgh'];
            mgh_corrected_file = [subj_dir,'/SurfaceProjections/equi_',hemisphere,'_whole_brain_',contrast,'_0p3_',dstr,'_CORRECTED.mgh'];
            mgh_file_copy = [subj_dir,'/SurfaceProjections/equi_',hemisphere,'_whole_brain_',contrast,'_0p3_',dstr,'.mgh'];
            
            % Use different paths for depth=0 (average across depths)
            if depth == 0
                mgh_file = [subj_dir,'/SurfaceProjections/',contrast,'_0p3_average_',hemisphere,'.mgh'];
                mgh_corrected_file = [subj_dir,'/SurfaceProjections/',contrast,'_0p3_average_',hemisphere,'_CORRECTED.mgh'];
                mgh_file_copy = [subj_dir,'/SurfaceProjections/',contrast,'_0p3_average_',hemisphere,'.mgh'];
            end
            
            % Handle file copying between locations if needed
            if exist(mgh_file) == 0 && exist(mgh_file_copy) == 1
                system(['cp ', mgh_file_copy, ' ', mgh_file]);
            end
            if exist(mgh_file) == 2 && exist(mgh_file_copy) == 0
                system(['cp ', mgh_file, ' ', mgh_file_copy]);
            end
            
            % Process file if it exists and hasn't been corrected yet
            if exist(mgh_file) == 2
                if exist(mgh_corrected_file) == 0
                    % Load MGH file
                    [vol, M, mr_parms, volsz] = load_mgh(mgh_file);
                    
                    % Identify and correct outliers
                    [outlidx, LTHRESH, UTHRESH, CENTER] = isoutlier(vol);
                    vol_corr = vol;
                    vol_corr(vol_corr < LTHRESH) = LTHRESH;  % Clip lower outliers
                    vol_corr(vol_corr > UTHRESH) = UTHRESH;  % Clip upper outliers
                    
                    % Save corrected volume
                    save_mgh(vol_corr, mgh_corrected_file, M, mr_parms);

                    % Generate visualization plots
                    surface_file = [subj_dir,'/surf/',hemisphere,'.pial'];
                    gg_newmat = gifti(surface_file);
                    
                    % Create plots for both uncorrected and corrected data
                    for type = 1:2
                        if type == 1
                            gg.cdata = single(vol);
                            string = 'uncorrected';
                        elseif type == 2
                            gg.cdata = single(vol_corr);
                            string = 'corrected';
                        end
                        
                        % Set up figure with consistent visualization parameters
                        rotation_angle = [90];
                        figure, plot(gg_newmat,gg);
                        view([rotation_angle 0])

                        % Configure lighting and colormap
                        hfig = light; 
                        set(hfig, 'position', [1 1 0.2]); 
                        lighting gouraud; 
                        material dull
                        caxis([display_range(1), display_range(end)]);
                        colormap(brewermap(256, '*Spectral'));
                        
                        % Add and configure colorbar
                        hfig = colorbar();
                        set(hfig,'XTick',display_range,'FontSize',10);
                        ylabel(hfig, contrast_names{c},'FontSize',11);
                        set(gcf,'color','w');
                        
                        % Save figure
                        figure_file_to_save_prefix = [plot_folder,'/Individual_figure_brewermap_',subj_dirs(s).name,'_',contrast,'_',hemisphere,'_',string];
                        figure_file_to_save = [figure_file_to_save_prefix,'.jpg'];
                        saveas(gcf,figure_file_to_save);
                    end
                end
            else
                display('mgh missing')
                mgh_file
            end
        end
        end

    end
end
