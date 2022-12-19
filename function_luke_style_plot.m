function function_luke_style_plot(templ_file, surface_file, valuelabels, valuetype, figure_file_to_save_prefix, projection_file_to_save, display_range)

    rotation_angles =  [90 -90]
    for ra = 1:length(rotation_angles)
        rotation_angle = rotation_angles(ra);

        %%% load template projection file
       % templ_file = [freesurfer_folder,'/Profiles/equi_',hemisphere,'_whole_brain_R1_0p3_0.5.mgh']; 
        [vol, M, mr_parms, volsz] = load_mgh(templ_file);
%         projection_file_to_save = [outdir_projections,'/Group_projection_',atlas,'_',contrasts{c},'_',hemisphere,'.mgh'];
%         if strcmp(hemisphere,'lh')
%             save_mgh(value_labels_lh, projection_file_to_save, M, mr_parms);
%             cdata = single(value_labels_lh);
%         else
%             save_mgh(value_labels_rh, projection_file_to_save, M, mr_parms);
%             cdata = single(value_labels_rh);
%         end
        save_mgh(valuelabels, projection_file_to_save, M, mr_parms);
        cdata = single(valuelabels);
        %% luke's code
        gg_newmat = gifti(surface_file);
        %cdata = single(mri.vol);
        gg.cdata = cdata';
        figure, plot(gg_newmat,gg);
        view([rotation_angle 0])
%             if strcmp(hemisphere,'lh')
%                 view([-90 0]); 
%             else
%                 view([90 0]); 
%             end
        hfig = light; set(hfig, 'position', [1 1 0.2]); lighting gouraud; material dull
        caxis([display_range(1), display_range(end)]);
        colormap(brewermap(256, '*Spectral'));
        hfig = colorbar();
        set(hfig,'XTick',display_range,'FontSize',10);
        ylabel(hfig, valuetype,'FontSize',11);
        set(gcf,'color','w');
        %figure_file_to_save = [resultdir,'/figures/Group_projection_',atlas,'_',contrasts{c},'_',hemisphere,'_',num2str(ra),'.jpg'];
        figure_file_to_save = [figure_file_to_save_prefix,'_',num2str(ra),'.jpg'];
        saveas(gcf,figure_file_to_save);
    end
    %%
end
