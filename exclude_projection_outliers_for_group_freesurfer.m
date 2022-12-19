close all
clear all
clc

addpath(genpath('freesurfer/6.0.0/ubuntu-xenial-amd64/matlab'))
addpath(genpath('Software/BrewerMap'))

analysis_folder = 'Cortical_analysis/group_freesurfer';
contrasts = {'R1','MTsat','R2s'};
contrast_names = {'R1 (s^{-1}) 7T','MT_{sat} (p.u.) 7T', 'R2* (s^{-1}) 7T''};

display_ranges = {[1.4 2],[3.5 4.5],[15 20]};
hemispheres = {'lh','rh'};

%%% make subject plot directory
plot_folder = [analysis_folder,'/individual_plots'];
system(['mkdir ',plot_folder]);

subj_dirs = dir([analysis_folder,'/Subj*']);
for s = 1:length(subj_dirs)
    subj_dir = [subj_dirs(s).folder,'/',subj_dirs(s).name];  

    for c = 1:length(contrasts)
      contrast = contrasts{c};
      display_range = display_ranges{c};
      for depth = [0 0.2 0.35 0.5 0.65 0.8]
          dstr = num2str(depth);
        for h = 1:2
            hemisphere = hemispheres{h};
           % mgh_file=[subj_dir,'/SurfaceProjections/',contrast,'_0p3_midcortical_',hemisphere,'.mgh'];
            mgh_file = [subj_dir,'/Profiles/equi_',hemisphere,'_whole_brain_',contrast,'_0p3_',dstr,'.mgh'];
            mgh_corrected_file = [subj_dir,'/SurfaceProjections/equi_',hemisphere,'_whole_brain_',contrast,'_0p3_',dstr,'_CORRECTED.mgh'];
            mgh_file_copy = [subj_dir,'/SurfaceProjections/equi_',hemisphere,'_whole_brain_',contrast,'_0p3_',dstr,'.mgh'];
            if depth == 0
                mgh_file = [subj_dir,'/SurfaceProjections/',contrast,'_0p3_average_',hemisphere,'.mgh'];
                mgh_corrected_file = [subj_dir,'/SurfaceProjections/',contrast,'_0p3_average_',hemisphere,'_CORRECTED.mgh'];
                mgh_file_copy = [subj_dir,'/SurfaceProjections/',contrast,'_0p3_average_',hemisphere,'.mgh'];
            end
            if exist(mgh_file) == 0 & exist(mgh_file_copy) == 1
                system(['cp ', mgh_file_copy, ' ', mgh_file]);
            end
            if exist(mgh_file) == 2 & exist(mgh_file_copy) == 0
                system(['cp ', mgh_file, ' ', mgh_file_copy]);
            end
            
            
            if exist(mgh_file) == 2
                if exist(mgh_corrected_file) == 0
                    %% manipulate
                    %mgh_file
                    %%% load file
                    [vol, M, mr_parms, volsz] = load_mgh(mgh_file);
                    %%% find outliers
                    %subplot(1,2,1)
                    %h = boxplot(vol);
                    [outlidx, LTHRESH, UTHRESH, CENTER] = isoutlier(vol);
                    vol_corr = vol;
                    %%% replace outliers with thresholds
                    vol_corr(vol_corr < LTHRESH) = LTHRESH;
                    vol_corr(vol_corr > UTHRESH) = UTHRESH;
                    %subplot(1,2,2)
                    %boxplot(vol_corr)
                    %%% save new file
                    save_mgh(vol_corr, mgh_corrected_file, M, mr_parms);

                    %% plots
    %                 
    %                 surface_file = [subj_dir,'/surf/',hemisphere,'.pial'];
    %                 gg_newmat = gifti(surface_file);
    %                 for type = 1:2
    %                     %subplot(h,2,type);
    %                     if type == 1
    %                         gg.cdata = single(vol);
    %                         string = 'uncorrected';
    %                     elseif type == 2
    %                         gg.cdata = single(vol_corr);
    %                         string = 'corrected';
    %                     end
    %                     rotation_angle = [90];
    % 
    %                     figure, plot(gg_newmat,gg);
    %                     view([rotation_angle 0])
    % 
    %                     hfig = light; set(hfig, 'position', [1 1 0.2]); lighting gouraud; material dull
    %                     caxis([display_range(1), display_range(end)]);
    %                     colormap(brewermap(256, '*Spectral'));
    %                     hfig = colorbar();
    %                     set(hfig,'XTick',display_range,'FontSize',10);
    %                     ylabel(hfig, contrast_names{c},'FontSize',11);
    %                     set(gcf,'color','w');
    %                     
    %                     figure_file_to_save_prefix = [plot_folder,'/Individual_figure_brewermap_',subj_dirs(s).name,'_',contrast,'_',hemisphere,'_',string];
    %                     figure_file_to_save = [figure_file_to_save_prefix,'.jpg'];
    %                     saveas(gcf,figure_file_to_save);
    %                 end
                end
            else
                display('mgh missing')
                mgh_file
            end
        end
        end

    end
end
