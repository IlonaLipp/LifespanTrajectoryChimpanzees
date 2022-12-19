function [globalprofile] = function_to_make_group_profiles(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, indir, outdir)
%%% finds all profile files that exist for a brain and averages them
   hemispheres = {'lh','rh'};
   types = {'','_partial'}
   for t = 1:2
     type = types{t};
     for con = 1:length(contrasts)
            contrast = contrasts{con};
            profcount = 0;
            for node = 1:nodenr
                
                for hem = 1:2
                      
                    hemisphere = hemispheres{hem};
                    filename_group_profile = [outdir,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.csv'];
                    filename_group_profile_raw = [outdir,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'_notnormalised.csv'];
                    filename_group_profile_plot = [outdir,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.png'];
                    filename_group_profile_plot_raw = [outdir,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'_notnormalised.png'];
                    filename_group_profile_troubleshooting_plot = [outdir,'/Group_profiles_troubleshooting_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.png'];
                    filename_group_medians = [outdir,'/Group_medians_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.csv'];

                    if exist(filename_group_profile_plot) == 0 || exist(filename_group_medians) == 0
                        %%% save profile as two column file with depth and
                        %%% value
                        filename_profiles = dir([indir,'/*_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.csv']);
                        for f = 1:length(filename_profiles)
                           data = csvread([filename_profiles(f).folder,'/', filename_profiles(f).name]);
                           %%% get chimp id
                           chimp_number(f) = str2double(filename_profiles(f).name(1:3));
                           for d = 1:length(depth_sampled)
                              idx = find(abs(data(:,1) - depth_sampled(d)) < .0001); %% equal does not work
                              all_profiles(d,f) = data(idx,2);
                           end
                           all_profiles_z(:,f) = (all_profiles(:,f) - mean(all_profiles(:,f)) ) ./ std(all_profiles(:,f)); 
                           if std(all_profiles(:,f)) == 0 %%% stupid workaround
                               all_profiles_z(:,f) = zeros(size(all_profiles_z(:,f)));
                           end
                        end

                        %%% average z-std profile, this makes sense to
                        %%% compare profiles trajectories across chimps
                        curr_profile = nanmedian(all_profiles_z'); %%% just in case
                        curr_profile_error = iqr(all_profiles_z');
                        csvwrite(filename_group_profile, [depth_sampled', curr_profile', curr_profile_error']);
                        
                        %%% not normalised, weighting higher the chimps
                        %%% with the highest values (so adults)
                        curr_profile_nn = nanmedian(all_profiles'); %%% just in case
                        curr_profile_nn_error = iqr(all_profiles');
                        csvwrite(filename_group_profile_raw, [depth_sampled', curr_profile_nn', curr_profile_nn_error']);
                        
                        
                        %% mean and skewness value per person (based on non z-standardised profiles)
                        regional_median = median(all_profiles);
                        regional_iqr = iqr(all_profiles);
                        %regional_skewness = skewness(all_profiles);
                        for b = 1:size(all_profiles,2)
                            %regional_skewness(1,b) = function_profile_skewness(depth_sampled, all_profiles(:,b));
                            %%% now new approach for calculating this as
                            %%% discussed with casey
                            regional_skewness(1,b) = function_profile_skewness(all_profiles(:,b));
                            %%% comparison parameter from paquola paper
                            regional_skewness_old_paquola(1,b) = skewness(all_profiles(:,b));
                            regional_skewness_my_old(1,b) = old_implementation_function_profile_skewness(depth_sampled, all_profiles(:,b));
                            if isnan(regional_skewness(1,b))
                                all_profiles(:,b)
                            end
                           legend_entry{b} = [num2str(chimp_number(b)), ': myold:', num2str(round(regional_skewness_my_old(1,b),2)) , '; CPold: ', num2str(round(regional_skewness_old_paquola(1,b),3)), '; CPnew: ', num2str(round(regional_skewness(1,b),3))];
                        end
                        
                        %% plot for sanity checking
                        colidx = find(~isnan(mean(all_profiles_z)));
                        %reduced_matrix = all_profiles_z(:,colidx);
                        reduced_matrix = all_profiles(:,colidx);
                        f = figure(300)
                            %set(f, 'visible','off');
                            set(f, 'position',[200 200 1500 500]);
                        %subplot(1,3,1)
                        subplot(1,2,1)
                            cmap = colormap(jet(length(colidx)));
                            %ax = axes('colororder',cmap);
                            plot(reduced_matrix)
                            ylabel(contrast)
                            %legend(num2str(chimp_number(colidx)'),'Location','NorthEastOutside')
                            legend(legend_entry(colidx),'Location','NorthEastOutside')
                        %subplot(1,3,2)  
                            %plot(nanmean(reduced_matrix'))
                            %ylabel(contrast)
                        %subplot(1,3,3)
                        subplot(1,2,2)
                            plot(nanmedian(reduced_matrix'),'b')
                            hold on
                            plot(nanmean(reduced_matrix'),'r')
                            ylabel(contrast)
                            legend({'median','mean'})
                        saveas(f,filename_group_profile_troubleshooting_plot);
                        close(f)
%                         %%% save in structure
%                         profcount = profcount + 1;
%                         all_node_profs(:,profcount) = curr_profile;
%                         
                        %%% plot
                        %filename_profile_plot = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.png'];
                        f = figure(1);
                            set(f, 'visible','off');
                            %group_skewness = function_profile_skewness(depth_sampled, curr_profile);
                            function_plot_profile(curr_profile, curr_profile_error, depth_sampled, contrast, node);
                            saveas(gcf,filename_group_profile_plot);
                            close(f)
                            
                        f = figure(1);
                            set(f, 'visible','off');
                            %group_skewness = function_profile_skewness(depth_sampled, curr_profile);
                            function_plot_profile(curr_profile_nn, curr_profile_nn_error, depth_sampled, contrast, node);
                            saveas(gcf,filename_group_profile_plot_raw);
                            close(f)    
                            

                        %figure()
                        %scatter(regional_skewness, regional_skewness_paquola)
                        %lsline
                        csvwrite(filename_group_medians, [chimp_number', regional_median', regional_iqr', regional_skewness', regional_skewness_old_paquola']);
                        
                    end
                end
            end

     end
        
     end
end
