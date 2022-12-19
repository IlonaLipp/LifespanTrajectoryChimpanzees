function function_to_plot_group_profiles_average_extended(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, outdir_group)

   hemispheres = {'lh','rh'};
   plotcols = {'r','b','g'};
   plotcols = {[123 15 15]./255, [20 15 123]./255, [19 92 45]./255};
   
    for node = 0:nodenr

        filename_profile_plot = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_ALL_2.png'];   
        if exist(filename_profile_plot) == 0
            f = figure(node+1);
            set(f,'Position',[100 100 1000 800]);
            %set(f,'units','normalized','Position',[0 0 .2 .45]);
            %set(f, 'visible','off');
            subpl_count = 0;
            for hem = 1:2
                for plottype = 1:2
                    subpl_count = subpl_count + 1;
                    for con = 2:3%length(contrasts)
                        contrast = contrasts{con};
                        hemisphere = hemispheres{hem};
                         if plottype == 1
                            filename_group_profile = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                         else
                            filename_group_profile = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'_partial.csv'];
                         end
                        profdata = csvread(filename_group_profile);

                        curr_profile = profdata(:,2); curr_profile_error = profdata(:,3); depth_sampled = profdata(:,1);

                        %%% plot 
                        subplot(2,2,subpl_count); 
                            hold on;
                            plot(curr_profile,depth_sampled(end:-1:1),'Color',plotcols{con},'LineWidth',2);
                            hold off
                            set(gca,'YTick',[0 max(depth_sampled)]);
                            set(gca,'YTickLabel',{'WM','pial'},'FontSize',12); %,'FontWeight','Bold');
                            xtickangle(45);
                            ylabel('cortical depth','FontSize',12) ; %,'FontWeight','Bold');
                            xlabel('z');
                            set(gcf,'color','white');
                            legend(contrasts(2:con),'Location','southwestoutside') %'southwest')
                            xlim([-2 2]); 
                            title(hemisphere)
                    end
                end
            end
            saveas(gcf,filename_profile_plot);
            close(f)
        end
    end
       
end
