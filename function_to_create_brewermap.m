function function_to_create_brewermap(surface_file, vol_to_plot, vol_name, title_name, figure_file_to_save_prefix, display_range)

    gg.cdata = single(vol_to_plot);            

    gg_newmat = gifti(surface_file); 

    rotation_angles =  [90 -90]
    for ra = 1:length(rotation_angles)

        rotation_angle = rotation_angles(ra);

        figure, plot(gg_newmat,gg);
        title(title_name);
        view([rotation_angle 0]);

        hfig = light; set(hfig, 'position', [1 1 0.2]); lighting gouraud; material dull

        %display_range = [prctile(vol_to_plot,10), prctile(vol_to_plot,90)]; %%% overwrrite from before;
        caxis([display_range(1), display_range(end)]);
        colormap(brewermap(256, '*Spectral'));
       % colormap()
        hfig = colorbar();
        set(hfig,'XTick',display_range,'FontSize',10);

        ylabel(hfig, vol_name,'FontSize',11);
        set(gcf,'color','w');

        figure_file_to_save = [figure_file_to_save_prefix,'_',num2str(ra),'.jpg'];
        try
            saveas(gcf,figure_file_to_save);
        catch
            display('fig cannot be saved')
        end
    end