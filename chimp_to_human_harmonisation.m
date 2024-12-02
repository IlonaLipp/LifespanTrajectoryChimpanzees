%% Chimpanzee to Human Brain Mapping Harmonisation Script
% This script performs harmonisation of cortical measurements between chimpanzee and human brain data.
% It analyzes R1, MTsat and R2s maps, comparing values across species and generating conversion factors.
%
% Key functions:
% 1. Loads and processes group-level data for both species
% 2. Performs regional analyses using BB38 atlas
% 3. Calculates conversion parameters between species
% 4. Applies conversion to individual chimpanzee data
%
% Dependencies:
% - FreeSurfer 6.0.0 MATLAB tools
% - BrewerMap for visualization
% - Custom functions: errorbarxy, function_wtlsc_line, function_to_create_brewermap

clear all
close all
clc

% Add required paths
addpath(genpath('freesurfer/6.0.0/ubuntu-xenial-amd64/matlab'))
rmpath(genpath('Primate_resources/HCPpipelines-master'))
addpath(genpath('Software/BrewerMap'))
addpath(fullfile(pwd, '../functions'));

resultdir = 'Cortical_analysis';

% Initialize variables
hemispheres = {'lh', 'rh', 'both'};
all_chimp_mean = [];
all_human_mean = [];
all_chimp_std = [];
all_human_std = [];
all_region_names = [];
figure()
mpm = 'R2s'

contrasts = {'R1','MTsat','R2s'};

% Load and summarize human demographic data
load_human_ages
mean(ages_from_database(1:15))
std(ages_from_database(1:15))
load_human_sex
females = sum(sex_from_database(1:15) == 1)
males = sum(sex_from_database(1:15) == 2)

adult_chimp_age = [ 17.5, 30, 34, 40, 43, 44, 45, 47, 52];
mean(adult_chimp_age)
std(adult_chimp_age)

for c = 1:length(contrasts)
    
    clear chimp_vals human_vals
    
    contrast = contrasts{c};
    residus = []; % Initialize residuals array
    both_hem_means_chimps = [];
    both_hem_means_humans = [];
    both_hem_std_chimps = [];
    both_hem_std_humans = [];
    
    reg_count = 0;
    
    % Process each hemisphere
    for h = 1:length(hemispheres)

        hemisphere = hemispheres{h};

        % Process individual hemispheres
        if strcmp(hemisphere, 'lh') || strcmp(hemisphere, 'rh')
            
            % Load chimp group average data
            chimpfile = ['Cortical_analysis/group_freesurfer/target_fsaverage/target_fsaverage_',hemisphere,'.main_just_adults_',contrast,'_0.5.mgh'];
            [volmat, M, mr_parms, volsz] = load_mgh(chimpfile);
            volmat_chimps = squeeze(volmat); % vertex x chimp matrix

            % Load human group average data
            if strcmp(contrast, 'R2s')
                humanfile = ['Cortical_analysis/group_freesurfer/humans/target_fsaverage/target_fsaverage_',hemisphere,'.main_R2s_WOLS_0.5.mgh'];
            else
                humanfile = ['Cortical_analysis/group_freesurfer/humans/target_fsaverage/target_fsaverage_',hemisphere,'.main_',contrast,'_0.5.mgh'];
            end
            [volmat, M, mr_parms, volsz] = load_mgh(humanfile);
            volmat_humans = squeeze(volmat); % vertex x human matrix

            % Load and process BB38 atlas information
            annotfile = ['Cortical_analysis/group_freesurfer/fsaverage/label/',hemisphere,'.BB38chimp.annot']
            [vertices, labeling, colortable] = read_annotation(annotfile);
            colors = colortable.table(1:end, 5);
            
            % Correct labeling indices
            clear labeling_corr
            for v = 1:length(labeling)
                if labeling(v) == 0
                    labeling_corr(v,1) = 0; 
                else
                    labeling_corr(v,1) = find(colors == labeling(v)) - 1; % -1 for BB38 
                end
            end

            % Extract regional values
            for r = 1:38
                vertidx = find(labeling_corr == r);
                reg_count = reg_count + 1;
                if length(vertidx) > 0
                   % Calculate median values per region for chimps
                   chimp_vals(:,reg_count) = median(volmat_chimps(vertidx,:),1)';
                    chimp_mean(r) = mean(chimp_vals(:,reg_count));
                    chimp_std(r) = std(chimp_vals(:,reg_count));
                   
                   % Calculate median values per region for humans
                   human_vals(:,reg_count) = median(volmat_humans(vertidx,:),1)';
                    human_mean(r) = mean(human_vals(:,reg_count));
                    human_std(r) = std(human_vals(:,reg_count));
                else
                   chimp_mean(r) = NaN;
                   chimp_std(r) = NaN;
                   human_mean(r) = NaN;
                   human_std(r) = NaN;
                end
            end
            
            % Stack hemisphere values
            both_hem_means_chimps = [both_hem_means_chimps, chimp_mean];
            both_hem_means_humans = [both_hem_means_humans, human_mean];
            both_hem_std_chimps = [both_hem_std_chimps, chimp_std];
            both_hem_std_humans = [both_hem_std_humans, human_std];
            chimp_means_to_use = chimp_mean;
            chimp_std_to_use = chimp_std;
            human_means_to_use = human_mean;
            human_std_to_use = human_std;
        else % For 'both' hemispheres case
            chimp_means_to_use = both_hem_means_chimps;
            chimp_std_to_use = both_hem_std_chimps;
            human_means_to_use = both_hem_means_humans;
            human_std_to_use = both_hem_std_humans;
            roi_names = [colortable.struct_names(2:end); colortable.struct_names(2:end)];
        end
        
        % Generate comparison plots
        figure()
            hold on
            scatter(chimp_means_to_use, human_means_to_use);
            xlabel(['chimp mean ', contrast]);
            ylabel(['human mean ', contrast]);
            
            % Add error bars
            xin = chimp_means_to_use;
            yin = human_means_to_use;
            uxin = chimp_std_to_use;
            uyin = human_std_to_use;
            errorbarxy(xin, yin, uxin, uyin);
            
            % Fit weighted total least squares line
            [a,b,alpha,p,chiopt,Cab,Calphap] = function_wtlsc_line(xin,yin,uxin,uyin);
            % a,b: usual straight line parameters; y=a*x+b
            refline(a,b);
            title([hemisphere, ': ', num2str(a) ' x chimp + ', num2str(b)]);
            
            % Calculate residuals
            residus = [residus; [yin - (a * xin + b)]'];
            
            % Save plot
            save_ylim = ylim;
            save_xlim = xlim;
            figure_file_to_save = [resultdir,'/figures/chimp_vs_human_',contrast,'_',hemisphere,'.jpg'];
            saveas(gcf,figure_file_to_save);

        % Generate surface maps for individual hemispheres
        if strcmp(hemisphere, 'lh') || strcmp(hemisphere, 'rh')
            % Calculate relative difference maps
            human_avg_map = mean(volmat_humans, 2); human_std_map = [std(volmat_humans')]';
            chimp_avg_map = mean(volmat_chimps, 2); chimp_std_map = [std(volmat_chimps')]';
            perc_change_chimps_to_humans = 100*(human_avg_map - chimp_avg_map)./(chimp_avg_map);
            surface_file = ['/data/pt_02101/results/Cortical_analysis/group_freesurfer/fsaverage/surf/',hemisphere,'.inflated'];
            display_range = [prctile(perc_change_chimps_to_humans,10), prctile(perc_change_chimps_to_humans, 90)];
            figure_file_to_save = [resultdir,'/figures/chimp_vs_human_',contrast,'_',hemisphere,'_relative_difference_hum_vs_chimps'];
            function_to_create_brewermap(surface_file, perc_change_chimps_to_humans, 'perc diff', contrast, figure_file_to_save, display_range);
        end
    end
    
    % Apply conversion to individual chimp data
    subj_dirs = dir(['Cortical_analysis/group_freesurfer/Subj*']);
    for s = 1:length(subj_dirs)
        subj_dir = [subj_dirs(s).folder,'/',subj_dirs(s).name];  
          for h = 1:2
                hemisphere = hemispheres{h};
                mgh_converted = [subj_dir,'/SurfaceProjections/',contrast,'_0p3_0.5_',hemisphere,'_converted_to_in_vivo_CORRECTED.mgh'];
                if exist(mgh_converted) == 0
                    mgh_corrected_file = [subj_dir,'/SurfaceProjections/equi_',hemisphere,'_whole_brain_',contrast,'_0p3_0.5_CORRECTED.mgh'];
                    if exist(mgh_corrected_file )
                        [volmat, M, mr_parms, volsz] = load_mgh(mgh_corrected_file);
                        volmat_corrected = a * volmat + b; % Apply conversion parameters
                    	save_mgh(volmat_corrected, mgh_converted, M, mr_parms);
                    end
                end
          end
    end
end
