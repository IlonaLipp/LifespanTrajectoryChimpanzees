function function_to_make_age_grouped_profiles(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, indir, outdir)

   load_chimp_ages
   %%% finds all profile files that exist for a brain and averages them
   hemispheres = {'lh','rh'};
     for con = 1:length(contrasts)
            contrast = contrasts{con};
            profcount = 0;
            for node = 1:nodenr
                for hem = 1:2
                      
                    hemisphere = hemispheres{hem};
                    filename_group_profile_plot = [outdir,'/Group_age_grouped_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.png'];

                    if 1 == 1 %exist(filename_group_profile_plot) == 0 
                        %%% save profile as two column file with depth and
                        %%% value
                        filename_profiles = dir([indir,'/*_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv']);
                        for f = 1:length(filename_profiles)
                           data = csvread([filename_profiles(f).folder,'/', filename_profiles(f).name]);
                           %%% get chimp id
                           chimp_number(f) = str2double(filename_profiles(f).name(1:3));
                           for d = 1:length(depth_sampled)
                              idx = find(abs(data(:,1) - depth_sampled(d)) < .0001); %% equal does not work
                              all_profiles(d,f) = data(idx,2);
                           end
                        end
                        
                        %%% 
                        ages_to_group = ages_from_database(chimp_number); 
                        
                        f = figure(1);
                        set(f, 'visible','off');
                        hold on;
                        colstouse = {'r','g','b'};
                        for group = 1:3
                            if group == 1
                                gridx = find(ages_to_group < 6);
                            elseif group == 2
                                gridx = find(ages_to_group >= 6 & ages_to_group < 17);
                            elseif group == 3
                                gridx = find(ages_to_group >= 17);
                            end
                            curr_profile = nanmedian(all_profiles(:,gridx)');
                            curr_profile_error = iqr(all_profiles(:,gridx)');
                            %se = shadedErrorBar(1:length(curr_profile), curr_profile, 0.5* curr_profile_error,colstouse{group},'transparent') %{'markerfacecolor','r'},0.25) 
                            %se.Annotation.LegendInformation.IconDisplayStyle = 'off';
                           % plot(curr_profile,'b-');
                            %plot(curr_profile,[colstouse{group},'-'],'LineWidth',3);
                            plot(depth_sampled',curr_profile',[colstouse{group},'-'],'LineWidth',3);
                        end
                        hold off;
                       % set(gca,'XTick',[1,length(depth_sampled)]);
                       % set(gca,'XTickLabel',{'pial','WM'},'FontSize',15,'FontWeight','Bold');
                        set(gca,'XTick',[0.05 0.25 0.5 0.75 max(depth_sampled)]);
                        set(gca,'XTickLabel',{'WM','25%','50%','75%','pial'},'FontSize',16); %,'FontWeight','Bold');
                        xtickangle(45);
                        xlabel('cortical depth','FontSize',15,'FontWeight','Bold');
                        ylabel(contrast);
                        title(['node ',sprintf('%.3d',node),' ',hemisphere]);
                        legend({'<6','6-17','17+'},'Location','NorthEastOutside');
                        set(gcf,'color','white');
                        saveas(gcf,filename_group_profile_plot);
                        close(f)    
                                 
                    end
                end
            end

     end
        
 end

