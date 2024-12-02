% This script generates group-level surface plots for FreeSurfer results comparing human and chimpanzee data
% It creates visualizations for different contrasts (R1, MTsat, R2s, thickness) across hemispheres
% Key features:
% - Processes both adult and child groups
% - Handles multiple contrasts and cortical depths 
% - Generates surface plots using BrewerMap colormaps
% - Saves results for both species comparisons

close all
clear all
clc

% Add required paths for FreeSurfer MATLAB tools and BrewerMap colormaps
addpath(genpath('freesurfer/6.0.0/ubuntu-xenial-amd64/matlab'))
addpath(genpath('Software/BrewerMap'))

% Set up analysis directories
analysis_folder = 'Cortical_analysis/group_freesurfer';
plot_folder = ['Cortical_analysis/group_freesurfer/Group_plots'];
system(['mkdir ', plot_folder]);

% Define analysis parameters
contrasts = {'R1','MTsat','R2s','thickness'};
contrast_names = {'R1 (s^{-1}) 7T','MT_{sat} (p.u.) 7T', 'R2* (s^{-1}) 7T', 'thickness (mm)'};

targets = {'fsaverage'} ; 
surfaces = {'inflated'}; % Can also use 'pial' or 'white'

for t = 1 
    target = targets{t};
    template_path = [analysis_folder,'/',target];
    hemispheres = {'lh','rh'};
   
    % Process different group effects
    effects = {'GroupAverageAdults','GroupAverageChildren','GroupAverage','Humans'};
    for e = [1 4] % Only process adults and humans
        % Set display ranges based on group type
        if strcmp(effects{e},'GroupAverageAdults') || strcmp(effects{e},'GroupAverageChildren')
             display_ranges = {[1.2 1.8],[3.7 5.3],[23.5 38.5],[2.2 3.3],[1.3 2],[17.5 28],[1.8 3.45]};
        elseif strcmp(effects{e},'GroupAverage')
             display_ranges = {[1.3 1.8],[3.6 5.1],[20.6 31.8],[2.4 3],[1.3 1.8],[17 22],[1.8 3.2]};
        else
             display_ranges = {[0 0], [0 0], [0 0], [0 0],[0 0], [0 0], [0 0]};
        end
        
        % Process each contrast type
        for c = 1:3 
          contrast = contrasts{c};

            % Handle different cortical depths
            for depth = [0 0.5] 
              if strcmp(contrast,'thickness')
                  depth = 0;
              else
                  depth = 0.5
              end

              dstr = num2str(depth);
              if depth == 0
                  dstr = 'avg'
              end
              
              % Process each hemisphere
              for h = 1:2
                hemisphere = hemispheres{h};
                
                % Load appropriate projection files based on effect type
                if strcmp(effects{e},'GroupAverageAdults')
                    % Load adult group files
                    projection_file_lh = [analysis_folder, '/target_',target,'/target_',target,'_lh.main_just_adults_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                    projection_file_rh = [analysis_folder, '/target_',target,'/target_',target,'_rh.main_just_adults_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                elseif strcmp(effects{e},'GroupAverageChildren')
                    % Load children group files
                    projection_file_lh = [analysis_folder, '/target_',target,'/target_',target,'_lh.main_just_children_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                    projection_file_rh = [analysis_folder, '/target_',target,'/target_',target,'_rh.main_just_children_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                elseif strcmp(effects{e},'GroupAverage')
                    % Load overall group average files
                    projection_file_lh = [analysis_folder, '/target_',target,'/target_',target,'_lh.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                    projection_file_rh = [analysis_folder, '/target_',target,'/target_',target,'_rh.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                elseif strcmp(effects{e},'Humans')
                    % Handle special case for R2s in humans
                    if strcmp(contrasts{c},'R2s')
                        projection_file_lh = [analysis_folder, '/humans/target_',target,'/target_',target,'_lh.main_R2s_WOLS_',dstr,'_sm03.osgm/beta.mgh'];
                        projection_file_rh = [analysis_folder, '/humans/target_',target,'/target_',target,'_rh.main_R2s_WOLS_',dstr,'_sm03.osgm/beta.mgh'];
                    else
                        projection_file_lh = [analysis_folder, '/humans/target_',target,'/target_',target,'_lh.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                        projection_file_rh = [analysis_folder, '/humans/target_',target,'/target_',target,'_rh.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                    end
                end                        

                % Load and process volume data
                [vol_lh, M, mr_parms, volsz] = load_mgh(projection_file_lh);
                [vol_rh, M, mr_parms, volsz] = load_mgh(projection_file_rh);
                vol_con = [vol_lh; vol_rh];
                volidx = vol_con~=0;

                % Calculate display ranges
                if strcmp(effects{e},'GroupAverage') || strcmp(effects{e},'GroupAverageAdults') ||strcmp(effects{e},'Humans')
                    display_range = [round(prctile(vol_con(volidx),2),1),round(prctile(vol_con(volidx),99), 2)];
                else
                    display_range = display_ranges{c};
                end

                % Load and save appropriate files for visualization
                if strcmp(effects{e},'GroupAverageAdults')
                    projection_file = [analysis_folder, '/target_',target,'/target_',target,'_',hemisphere,'.main_just_adults_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                    new_fname =  ['Downloads/Supplementary_Chimpanzees_map_',contrast,'_',hemisphere,'.mgh']
                    system(['cp ', projection_file, ' ', new_fname])
                elseif  strcmp(effects{e},'GroupAverageChildren')
                   projection_file = [analysis_folder, '/target_',target,'/target_',target,'_',hemisphere,'.main_just_children_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                elseif strcmp(effects{e},'GroupAverage')
                   projection_file = [analysis_folder, '/target_',target,'/target_',target,'_',hemisphere,'.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                elseif strcmp(effects{e},'Humans')
                     if strcmp(contrasts{c},'R2s')
                        projection_file = [analysis_folder, '/humans/target_',target,'/target_',target,'_',hemisphere,'.main_R2s_WOLS_',dstr,'_sm03.osgm/beta.mgh'];
                     else
                        projection_file = [analysis_folder, '/humans/target_',target,'/target_',target,'_',hemisphere,'.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                     end
                     new_fname =  ['Downloads/Supplementary_Human_map_',contrast,'_',hemisphere,'.mgh']
                    system(['cp ', projection_file, ' ', new_fname])
                end
                [vol, M, mr_parms, volsz] = load_mgh(projection_file);
                    
                % Generate visualization if volume exists   
                if exist('vol') == 1
                        for su = 1
                            % Set up plotting parameters
                            vol_to_plot = vol;
                            if strcmp(effects{e},'AgeEffect')
                                vol_to_plot(pval > .05) = 0;
                            end
                        
                            surface = surfaces{su} 
                            surface_file = [template_path,'/surf/',hemisphere,'.',surface];
    
                            % Set volume name based on effect type
                            if strcmp(effects{e},'AgeEffect')
                                vol_name =  ['Age model slope ', contrast_names{c}];
                            else
                                vol_name = contrast_names{c};
                            end

                            % Generate output filenames
                            figure_file_to_save_prefix = [plot_folder,'/Group_figure_brewermap_',target,'_',effects{e},'_',contrast,'_',hemisphere,'_',dstr,'_',surface];
                            title_name = [effects{e},' ',dstr];
                            
                            % Create and save visualization
                            function_to_create_brewermap(surface_file, vol_to_plot, vol_name, title_name, figure_file_to_save_prefix, display_range);
                        end
                end
              end
        end
    end
end
