clear all
close all
clc

addpath(genpath('freesurfer/6.0.0/ubuntu-xenial-amd64/matlab'))
addpath(genpath('Software/BrewerMap'))

contrasts = {'R1','MTsat','R2s','thickness'};
hemispheres = {'lh','rh'};

for c = 1:length(contrasts)
    
    contrast = contrasts{c};
              
    for hem = 1:length(hemispheres)
        
        hemisphere = hemispheres{hem};

        if strcmp(contrast, 'thickness')     
            cd = 'avg';
        else
            cd = '0.5';
        end
        
         foldername = ['/data/pt_02101/results/Cortical_analysis/group_freesurfer/chimps_vs_humans/target_fsaverage/target_fsaverage_',hemisphere,'.',contrast,'_',cd,'_sm03.osgm/humans_vs_chimps'];

        surface_file = ['/data/pt_02101/results/Cortical_analysis/group_freesurfer/chimps_vs_humans/fsaverage/surf/', hemisphere,'.inflated'];
     
        %%% rename and copy cluster stats name
        orig_fname =  [foldername, '/cache.th30.abs.sig.cluster.summary'];
        new_fname =  ['Downloads/cluster_stats_',contrast,'_',hemisphere,'_',cd]
        system(['cp ', orig_fname, ' ', new_fname])
        
        %%% rename and copy mgh filename
        orig_fname =  [foldername, '/cache.th30.abs.sig.cluster.mgh'];
        new_fname =  ['Downloads/Supplementary_Humans_vs_chimpanzees_clusters_',contrast,'_',hemisphere,'.mgh']
        system(['cp ', orig_fname, ' ', new_fname])

        vars_to_plot = {'z', 'clusters', 'z_in_clusters', 'chimps_converted_only', 'humans_only'};
        for v = 1:length(vars_to_plot) %%% badly written piece of code
            var_to_plot = vars_to_plot{v};
           if strcmp(var_to_plot, 'z')
               [vol_to_plot, M, mr_parms, volsz] = load_mgh([foldername, '/z.mgh']);
               display_range = [-4 4];
           elseif strcmp(var_to_plot, 'clusters')
               [vol_to_plot, M, mr_parms, volsz] = load_mgh([foldername, '/cache.th30.abs.sig.cluster.mgh']);
               display_range = [-4 4];
           elseif strcmp(var_to_plot, 'z_in_clusters')
               [vol_to_plot, M, mr_parms, volsz] = load_mgh([foldername, '/z.mgh']);
               [clu, M, mr_parms, volsz] = load_mgh([foldername, '/cache.th30.abs.sig.cluster.mgh']);
               vol_to_plot(clu == 0) = 0; %%% threshold
               display_range = [min(vol_to_plot) max(vol_to_plot)];
               display_range = [-4 4];
           elseif strcmp(var_to_plot, 'chimps_converted_only')
               [vol_to_plot, M, mr_parms, volsz] = load_mgh(['Cortical_analysis/group_freesurfer/chimps_vs_humans/target_fsaverage/target_fsaverage_',hemisphere,'.',contrast,'_',cd,'_sm03.osgm/only_chimps/gamma.mgh']);
               display_range = [prctile(vol_to_plot,10) prctile(vol_to_plot,90)];
           elseif strcmp(var_to_plot, 'humans_only')
               [vol_to_plot, M, mr_parms, volsz] = load_mgh(['Cortical_analysis/group_freesurfer/chimps_vs_humans/target_fsaverage/target_fsaverage_',hemisphere,'.',contrast,'_',cd,'_sm03.osgm/only_humans/gamma.mgh']);
               %%% do not reset display range now
           end
           title_name = var_to_plot;
           vol_name = [contrast, ' ', var_to_plot];
           figure_file_to_save_prefix = ['Cortical_analysis/group_freesurfer/Group_plots/Group_figure_brewermap_fsaverage_humans_vs_chimps_',contrast,'_',hemisphere,'_',cd,'_', var_to_plot];
           function_to_create_brewermap(surface_file, vol_to_plot, vol_name, title_name, figure_file_to_save_prefix, display_range);
           
        end
            
    end
end
