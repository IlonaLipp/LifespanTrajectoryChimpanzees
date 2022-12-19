clear all
close all
clc

addpath(genpath('freesurfer/6.0.0/ubuntu-xenial-amd64/matlab'))
rmpath(genpath('Primate_resources/HCPpipelines-master'))
addpath(genpath('Software/BrewerMap'))

resultdir = 'Cortical_analysis';

hemispheres = {'lh', 'rh', 'both'};
all_chimp_mean = [];
all_human_mean = [];
all_chimp_std = [];
all_human_std = [];
all_region_names = [];
figure()
mpm = 'R2s'

contrasts = {'R1','MTsat','R2s'};

%%% put together stats for humans
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
    residus = []; %%% initiate
    both_hem_means_chimps = [];
    both_hem_means_humans = [];
    both_hem_std_chimps = [];
    both_hem_std_humans = [];
    
    reg_count = 0;
    
    for h = 1:length(hemispheres)

        hemisphere = hemispheres{h};

        if strcmp(hemisphere, 'lh') || strcmp(hemisphere, 'rh')
            
            %% chimps group average
            %%% group MTsat data for adult chimps
            chimpfile = ['Cortical_analysis/group_freesurfer/target_fsaverage/target_fsaverage_',hemisphere,'.main_just_adults_',contrast,'_0.5.mgh'];
            [volmat, M, mr_parms, volsz] = load_mgh(chimpfile);
            volmat_chimps = squeeze(volmat); %%% should be vertex x chimp

            %% humans group average
            %%% group MTsat data for adult humans
            if strcmp(contrast, 'R2s')
                humanfile = ['Cortical_analysis/group_freesurfer/humans/target_fsaverage/target_fsaverage_',hemisphere,'.main_R2s_WOLS_0.5.mgh'];
            else
                humanfile = ['Cortical_analysis/group_freesurfer/humans/target_fsaverage/target_fsaverage_',hemisphere,'.main_',contrast,'_0.5.mgh'];
            end
            [volmat, M, mr_parms, volsz] = load_mgh(humanfile);
            volmat_humans = squeeze(volmat); %%% should be vertex x chimp

            %% get chimp atlas information
            annotfile = ['Cortical_analysis/group_freesurfer/fsaverage/label/',hemisphere,'.BB38chimp.annot']
            [vertices, labeling, colortable] = read_annotation(annotfile);
            colors = colortable.table(1:end, 5);
            %%% correct labelling
            clear labeling_corr
            for v = 1:length(labeling)
                if labeling(v) == 0
                    labeling_corr(v,1) = 0; 
                else
                    labeling_corr(v,1) = find(colors == labeling(v)) - 1; %%% -1 for bb38 
                end
            end

            %% extract values
            %%% we now have a labeling_corr vector with the region idx
            %%% get chimp mean values that correspond to luke values
            for r = 1:38 %length(regionnames)
                %roi = regionnames{r};
                %roiidx = find(strcmp(colortable.struct_names,roi));
                %vertidx = find(labeling_corr == roiidx);
                vertidx = find(labeling_corr == r);
                reg_count = reg_count + 1;
                if length(vertidx) > 0
                   %chimp_mean(r) = mean(vol(vertidx));
                   chimp_vals(:,reg_count) = median(volmat_chimps(vertidx,:),1)'; %%% 1 value per chimp
                    chimp_mean(r) = mean(chimp_vals(:,reg_count));
                    chimp_std(r) = std(chimp_vals(:,reg_count));
   %                chimp_mean(r) = median(chimp_vals);
   %               chimp_std(r) = iqr(chimp_vals); 
                   %%%
                   human_vals(:,reg_count) = median(volmat_humans(vertidx,:),1)'; %%% 1 value per chimp
                    human_mean(r) = mean(human_vals(:,reg_count));
                    human_std(r) = std(human_vals(:,reg_count));
   %                human_mean(r) = median(human_vals);
   %                human_std(r) = iqr(human_vals);     
                else
                   chimp_mean(r) = NaN;
                   chimp_std(r) = NaN;
                   human_mean(r) = NaN;
                   human_std(r) = NaN;
                end
            end
            %%% stack values of lh and rh and set "used" to original values
            both_hem_means_chimps = [both_hem_means_chimps, chimp_mean];
            both_hem_means_humans = [both_hem_means_humans, human_mean];
            both_hem_std_chimps = [both_hem_std_chimps, chimp_std];
            both_hem_std_humans = [both_hem_std_humans, human_std];
            chimp_means_to_use = chimp_mean;
            chimp_std_to_use = chimp_std;
            human_means_to_use = human_mean;
            human_std_to_use = human_std;
        else %%% set used to stacked values
            chimp_means_to_use = both_hem_means_chimps;
            chimp_std_to_use = both_hem_std_chimps;
            human_means_to_use = both_hem_means_humans;
            human_std_to_use = both_hem_std_humans;
            roi_names = [colortable.struct_names(2:end); colortable.struct_names(2:end)];
        end
        

        %% comparison 
        %%% scatter
        figure()
            hold on
            %idx = find(chimp_mean > 2); 
            %idx = 1:length(chimp_means_to_use)
            scatter(chimp_means_to_use, human_means_to_use);
            xlabel(['chimp mean ', contrast]);
            ylabel(['human mean ', contrast]);
            
            xin = chimp_means_to_use; %(idx);
            yin = human_means_to_use; %(idx);
            uxin = chimp_std_to_use; %(idx);
            uyin = human_std_to_use; %(idx);
            errorbarxy(xin, yin, uxin, uyin);
            %%% implement this special fit
            [a,b,alpha,p,chiopt,Cab,Calphap] = function_wtlsc_line(xin,yin,uxin,uyin);
            % a,b: usual straight line parameters; y=a*x+b
            refline(a,b);
            title([hemisphere, ': ', num2str(a) ' x chimp + ', num2str(b)]);
            %%% calculate confidence interval of slope
            %%% calculate residuals
            residus = [residus; [yin - (a * xin + b)]']; %%% create vector with residuasl for lh; rh
            %%% save plot
            save_ylim = ylim;
            save_xlim = xlim;
            figure_file_to_save = [resultdir,'/figures/chimp_vs_human_',contrast,'_',hemisphere,'.jpg'];
            saveas(gcf,figure_file_to_save);


        if strcmp(hemisphere, 'lh') || strcmp(hemisphere, 'rh')
            %% make whole brain relative difference maps, but how to decide what is significant?
            human_avg_map = mean(volmat_humans, 2); human_std_map = [std(volmat_humans')]';
            chimp_avg_map = mean(volmat_chimps, 2); chimp_std_map = [std(volmat_chimps')]';
            perc_change_chimps_to_humans = 100*(human_avg_map - chimp_avg_map)./(chimp_avg_map);
            surface_file = ['/data/pt_02101/results/Cortical_analysis/group_freesurfer/fsaverage/surf/',hemisphere,'.inflated'];
            display_range = [prctile(perc_change_chimps_to_humans,10), prctile(perc_change_chimps_to_humans, 90)];
            figure_file_to_save = [resultdir,'/figures/chimp_vs_human_',contrast,'_',hemisphere,'_relative_difference_hum_vs_chimps'];
            function_to_create_brewermap(surface_file, perc_change_chimps_to_humans, 'perc diff', contrast, figure_file_to_save, display_range);
        end
    end
    
    %% once this is done we have the slope for conversion, done for "both"
    %%% now i need to read and write out the individual subjects mgh files
    %%% just for chimps
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
                        volmat_corrected = a * volmat + b; %%% these are saved from the loop run for "both"
                    	save_mgh(volmat_corrected, mgh_converted, M, mr_parms);
                    end
                end
          end
    end
end
