 
function function_to_get_depth_vals(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, outdir, brain_undersc_id)
 
     hemispheres = {'lh','rh'};
     for con = 1:length(contrasts)
            contrast = contrasts{con};

            for node = 1:nodenr
                
                for hem = 1:2
                      
                    hemisphere = hemispheres{hem};

                    filename_profile = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                    filename_profile_plot = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.png'];

                    clear curr_profile curr_profile_error
                    
                    if exist(filename_profile_plot) == 0
                    
                        for d = 1:length(depth_sampled)

                                depth = depth_sampled(d);                 

                                projectionfile = [freesurfer_folder,'/Profiles/equi_',hemisphere,'_',atlas,'-',sprintf('%.3d',node),'_',contrast,'_0p3_',sprintf('%.02f',depth),'.mgh']
                                if exist(projectionfile) == 0
                                      projectionfile = [freesurfer_folder,'/Profiles/equi_',hemisphere,'_',atlas,'-',sprintf('%.3d',node),'_',contrast,'_0p3_',sprintf('%.01f',depth),'.mgh'];   
                                end 
                                if exist(projectionfile) == 2
                                    vals = load_mgh(projectionfile);
                                    curr_profile(d) = nanmedian(vals(vals~=0));
                                    curr_profile_error(d) = iqr(vals(vals~=0));
                                else %%% if no contrast file availale
                                    curr_profile(d) = NaN;
                                    curr_profile_error(d) = NaN;
                                end
                                %%% exclude R1 for a few
                                if strcmp(contrast(1:2),'R1')
                                    if strcmp(brain_undersc_id(1:3), '004') || strcmp(brain_undersc_id(1:3), '032') || strcmp(brain_undersc_id(1:3), '033')
                                        curr_profile(d) = NaN;
                                        curr_profile_error(d) = NaN;
                                    end
                                end
                                
                        end

                        %%% save profile as two column file with depth and
                        %%% value
                        csvwrite(filename_profile, [depth_sampled', curr_profile', curr_profile_error']);
                        %%% plot
                        f = figure(1)
                            set(f, 'visible','off');
                            function_plot_profile(curr_profile, curr_profile_error, depth_sampled, contrast, node)
                            saveas(gcf,filename_profile_plot);
                            close(f)
                    end
                end
             end
     end
     
     %%% global profiles
     hemispheres = {'lh','rh'};
     for con = 1:length(contrasts)
            contrast = contrasts{con};
            filename_globprofile_plot = [outdir,'/',brain_undersc_id,'_',atlas,'_global_',contrast,'.png'];
            filename_globprofile = [outdir,'/',brain_undersc_id,'_',atlas,'_global_',contrast,'.csv'];
            contrast = contrasts{con};
            all_profiles.(contrast) = [];
            for node = 1:nodenr
                for hem = 1:2
                    hemisphere = hemispheres{hem};
                    filename_profile = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                    profdata = csvread(filename_profile);
                    curr_profile = profdata(:,2);
                    curr_profile_z = (curr_profile - mean(curr_profile)) ./ std(curr_profile);
                    all_profiles.(contrast) = [all_profiles.(contrast), curr_profile_z];
                end
            end
           % global_profile.(contrast) = nanmean(all_profiles.(contrast), 2); %%% global profile is average of all
           % global_profile_error.(contrast) = nanstd(all_profiles.(contrast)')';
            global_profile.(contrast) = nanmedian(all_profiles.(contrast), 2); %%% global profile is average of all
            global_profile_error.(contrast) = iqr(all_profiles.(contrast)')';
%             if con == 3 & isnan(global_profile.(contrast)(1))
%                1 + 1
%                
%             end
            if exist(filename_globprofile_plot) == 0
                %%% plot global profile
                f = figure(1);
                    set(f, 'visible','off');
                    function_plot_profile(global_profile.(contrast), global_profile_error.(contrast), depth_sampled, contrast, 99)
                    saveas(gcf,filename_globprofile_plot);
                    close(f)
            end
            csvwrite(filename_globprofile, [depth_sampled', global_profile.(contrast), global_profile_error.(contrast)]);
     end
     %%% partial profiles: now regress out global profile
     hemispheres = {'lh','rh'};
     for con = 1:length(contrasts)
            contrast = contrasts{con};
            all_profiles.(contrast) = [];
            for node = 1:nodenr
                for hem = 1:2
                    hemisphere = hemispheres{hem};
                    filename_profile = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                    filename_partprofile = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'_partial.csv'];
                    filename_partprofile_plot = [outdir,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'_partial.png'];
                    if exist(filename_partprofile) == 0 || exist(filename_partprofile_plot) == 0
                        profdata = csvread(filename_profile);
                        curr_profile = profdata(:,2); 
                        curr_profile_error = profdata(:,3); %depth_sampled = profdata(:,1);
                        curr_profile_z = (curr_profile - mean(curr_profile)) ./ std(curr_profile);
                        [b bint r] = regress(curr_profile_z,[global_profile.(contrast),ones(length(curr_profile),1)]);
                       % [b bint r] = regress(curr_profile_z,[depth_sampled',ones(length(curr_profile),1)]);
                        partial_profile = r;
                        csvwrite(filename_partprofile, [depth_sampled', partial_profile, curr_profile_error]);
                        %%% plot
                        f = figure(1)
                            set(f, 'visible','off');
                            function_plot_profile(partial_profile, 0*curr_profile_error, depth_sampled, contrast, node)
                            saveas(gcf,filename_partprofile_plot);
                            close(f)
                    end
                end
            end
     end
end