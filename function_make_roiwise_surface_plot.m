function function_make_roiwise_surface_plot(variable_vec, variable_name, plot_perc, projection_file_to_save, figure_file_to_save_prefix)

    freesurfer_folder = 'Cortical_analysis/group_freesurfer/fsaverage';
    atlas_folder = 'LS_FS_atlas/BB38chimp_forSharing_Ilona';
    
    %volidx = find(variable_vec ~= 0);
    volidx = 1:length(variable_vec);
    display_range = [round(prctile(variable_vec(volidx),plot_perc(1)),2),round(prctile(variable_vec(volidx),plot_perc(2)), 2)];
    
    hemispheres = {'lh','rh'};
    for hem = 1:2
        
        hemisphere = hemispheres{hem}
        if hem == 1
            values_to_assign = variable_vec(1:38);
        else
            values_to_assign = variable_vec(39:76);
        end

        [vertices, labeling, colortable] = read_annotation([freesurfer_folder '/label/',hemisphere,'.BB38chimp.annot']);
        surface_file = ([freesurfer_folder '/surf/',hemisphere,'.white']);
        colors = colortable.table(1:end, 5);
           %%% correct labelling
            clear labeling_corr
            for v = 1:length(labeling)
                    labeling_corr(v,1) = find(colors == labeling(v))-1; 
            end

           %%% assign values
           clear val_vec
            for label = 1:length(labeling_corr)
                if labeling_corr(label) > 0
                    val_vec(label) = values_to_assign(labeling_corr(label));
                else
                    val_vec(label) = 0;
                end 
            end

        templ_file = [freesurfer_folder '/surf/',hemisphere,'.white.avg.area.mgh']
        [vol, M, mr_parms, volsz] = load_mgh(templ_file);
        try
            save_mgh(val_vec, projection_file_to_save, M, mr_parms);
        end
%         rotation_angles =  [90 -90]
% 
%         gg.cdata = single(val_vec + 2*rand(length(val_vec)));
          gg.cdata = single(val_vec)'; % + 0.01*rand(length(val_vec)));
         gg_newmat = gifti(surface_file); 
         
         %surf_conv = convert_surface(surface_file);
         % h = plot_hemispheres(gg.cdata,surf_conv) ;

          rotation_angles = [90 -90];
        for ra = 1:length(rotation_angles)
            rotation_angle = rotation_angles(ra);

            colormap(brewermap(256, '*Spectral'));
            figure, plot(gg_newmat,gg); 
            try
                caxis([display_range(1), display_range(end)]);
            catch
                display_range = [min(val_vec), max(val_vec)];
                caxis([display_range(1), display_range(end)]);
            end
            %title([effects{e},' ',dstr])
            view([rotation_angle 0])

            hfig = light; set(hfig, 'position', [1 1 0.2]); lighting gouraud; material dull
            colormap(brewermap(256, '*Spectral'));
            hfig = colorbar();
            
            try
            set(hfig,'XTick',display_range,'FontSize',10);
            end
            ylabel(hfig,variable_name,'FontSize',11,'interpreter','none');
            set(gcf,'color','w');
            
          
            figure_file_to_save = [figure_file_to_save_prefix,'_',num2str(ra),'_',hemisphere,'.jpg'];
            try %%% just so that it produces the figures even with incorrect path
                saveas(gcf,figure_file_to_save);
            end
        end
    end
end
