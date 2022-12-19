close all
clear all
clc


addpath(genpath('freesurfer/6.0.0/ubuntu-xenial-amd64/matlab'))
addpath(genpath('Software/BrewerMap'))

analysis_folder = 'Cortical_analysis/group_freesurfer';
plot_folder = ['Cortical_analysis/group_freesurfer/Group_plots'];
system(['mkdir ', plot_folder]);

contrasts = {'R1','MTsat','R2s','thickness'};
contrast_names = {'R1 (s^{-1}) 7T','MT_{sat} (p.u.) 7T', 'R2* (s^{-1}) 7T', 'thickness (mm)'};

targets = {'fsaverage'} ; 
surfaces = {'pial','inflated','white'};
surfaces = {'inflated'};

for t = 1 
    
    target = targets{t};
    template_path = [analysis_folder,'/',target];

    hemispheres = {'lh','rh'};
   
    effects = {'GroupAverageAdults','GroupAverageChildren','GroupAverage','Humans'};
    for e = [1 4] %1:length(effects)
        if strcmp(effects{e},'GroupAverageAdults') || strcmp(effects{e},'GroupAverageChildren')
             display_ranges = {[1.2 1.8],[3.7 5.3],[23.5 38.5],[2.2 3.3],[1.3 2],[17.5 28],[1.8 3.45]};
        elseif strcmp(effects{e},'GroupAverage')
             display_ranges = {[1.3 1.8],[3.6 5.1],[20.6 31.8],[2.4 3],[1.3 1.8],[17 22],[1.8 3.2]};
        else
             display_ranges = {[0 0], [0 0], [0 0], [0 0],[0 0], [0 0], [0 0]};
        end
        
        for c = 1:3 
          contrast = contrasts{c};

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
              
              for h = 1:2

                hemisphere = hemispheres{h};
                
                if strcmp(effects{e},'GroupAverageAdults')
                    %%% calculate optimal display ranges (same for left and
                    %%% right)
                    projection_file_lh = [analysis_folder, '/target_',target,'/target_',target,'_lh.main_just_adults_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                    projection_file_rh = [analysis_folder, '/target_',target,'/target_',target,'_rh.main_just_adults_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                elseif strcmp(effects{e},'GroupAverageChildren')
                    %%% calculate optimal display ranges (same for left and
                    %%% right)
                    projection_file_lh = [analysis_folder, '/target_',target,'/target_',target,'_lh.main_just_children_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                    projection_file_rh = [analysis_folder, '/target_',target,'/target_',target,'_rh.main_just_children_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];

                elseif strcmp(effects{e},'GroupAverage')
                    %%% calculate optimal display ranges (same for left and
                    %%% right)
                    projection_file_lh = [analysis_folder, '/target_',target,'/target_',target,'_lh.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                    projection_file_rh = [analysis_folder, '/target_',target,'/target_',target,'_rh.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                elseif strcmp(effects{e},'Humans')
                    if strcmp(contrasts{c},'R2s') %%% they are called differently
                        projection_file_lh = [analysis_folder, '/humans/target_',target,'/target_',target,'_lh.main_R2s_WOLS_',dstr,'_sm03.osgm/beta.mgh'];
                        projection_file_rh = [analysis_folder, '/humans/target_',target,'/target_',target,'_rh.main_R2s_WOLS_',dstr,'_sm03.osgm/beta.mgh'];
                    else
                        projection_file_lh = [analysis_folder, '/humans/target_',target,'/target_',target,'_lh.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                        projection_file_rh = [analysis_folder, '/humans/target_',target,'/target_',target,'_rh.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                    end
                end                        
                [vol_lh, M, mr_parms, volsz] = load_mgh(projection_file_lh);
                [vol_rh, M, mr_parms, volsz] = load_mgh(projection_file_rh);
                vol_con = [vol_lh; vol_rh];
                volidx = vol_con~=0;

                %display_range = [round(10*prctile(vol_con(volidx),10))/10,round(10*prctile(vol_con(volidx),90))/10];
                %%% set display range to 1st and 99th percentile for
                %%% average maps to start with
                if strcmp(effects{e},'GroupAverage') || strcmp(effects{e},'GroupAverageAdults') ||strcmp(effects{e},'Humans')
                    display_range = [round(prctile(vol_con(volidx),2),1),round(prctile(vol_con(volidx),99), 2)];
                else
                    display_range = display_ranges{c};
                end

                %%% looad correct for now
                if strcmp(effects{e},'GroupAverageAdults')
                    projection_file = [analysis_folder, '/target_',target,'/target_',target,'_',hemisphere,'.main_just_adults_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                    new_fname =  ['Downloads/Supplementary_Chimpanzees_map_',contrast,'_',hemisphere,'.mgh']
                    system(['cp ', projection_file, ' ', new_fname])
                elseif  strcmp(effects{e},'GroupAverageChildren')
                   projection_file = [analysis_folder, '/target_',target,'/target_',target,'_',hemisphere,'.main_just_children_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                elseif strcmp(effects{e},'GroupAverage')
                   projection_file = [analysis_folder, '/target_',target,'/target_',target,'_',hemisphere,'.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                elseif strcmp(effects{e},'Humans')
                     if strcmp(contrasts{c},'R2s') %%% they are called differently
                        projection_file = [analysis_folder, '/humans/target_',target,'/target_',target,'_',hemisphere,'.main_R2s_WOLS_',dstr,'_sm03.osgm/beta.mgh'];
                     else
                        projection_file = [analysis_folder, '/humans/target_',target,'/target_',target,'_',hemisphere,'.main_',contrast,'_',dstr,'_sm03.osgm/beta.mgh'];
                     end
                     new_fname =  ['Downloads/Supplementary_Human_map_',contrast,'_',hemisphere,'.mgh']
                    system(['cp ', projection_file, ' ', new_fname])
                end
                [vol, M, mr_parms, volsz] = load_mgh(projection_file);
                    
                   
                if exist('vol') == 1

                        for su = 1
                            
                            %%% set parameters for plotting function
                            
                            vol_to_plot = vol;
                            if strcmp(effects{e},'AgeEffect')
                                vol_to_plot(pval > .05) = 0;
                            end
                        
                            surface = surfaces{su} 
                            surface_file = [template_path,'/surf/',hemisphere,'.',surface];
    
                            if strcmp(effects{e},'AgeEffect')
                                vol_name =  ['Age model slope ', contrast_names{c}];
                            else
                                vol_name = contrast_names{c};
                            end
                            figure_file_to_save_prefix = [plot_folder,'/Group_figure_brewermap_',target,'_',effects{e},'_',contrast,'_',hemisphere,'_',dstr,'_',surface];
                            title_name = [effects{e},' ',dstr];
                            
                            %%% function
                            function_to_create_brewermap(surface_file, vol_to_plot, vol_name, title_name, figure_file_to_save_prefix, display_range);
                            

                        end
                end
              end
        end
    end
end
