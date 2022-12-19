 
function function_to_get_depth_vals_test2(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, outdir, brain_undersc_id)
 
     hemispheres = {'lh','rh'};
     
     
     for con = 1:length(contrasts)
            
           contrast = contrasts{con};
            
            for hem = 1:2
                    
                                    
                    clear all_vals new_vals
                    hemisphere = hemispheres{hem};
                    
                    %%% check if already done
                    donefiles = dir([outdir,'/TEST',brain_undersc_id,'_',atlas,'_*_',hemisphere,'_',contrast,'.csv']);
                    
                    if length(donefiles) < 38

                        %% get label
                        old_annotation_file = [freesurfer_folder,'/label/',hemisphere,'.BB38chimp.annot'];     
                        [vertices, label, colortable] = read_annotation(old_annotation_file);
                        BB38chimpIL = read_my_new_ctab('/data/tu_lippi/Software/Primate_resources/LS_FS_atlas/BB38chimp_forSharing_Ilona/BB38chimpIL.annot.ctab'); %%% script to open ctab file i created
                        roicodes = colortable.table(2:end,5); %%% code for each roi

                        %% get thickness
                        thickness_file = [freesurfer_folder,'/surf/',hemisphere,'.thickness'];
                        thicknessvals = read_curv(thickness_file);

                        %% get profile for each vertex
                        for d = 1:length(depth_sampled)
                            depth = depth_sampled(d);                 
                            projectionfile = [freesurfer_folder,'/Profiles/equi_',hemisphere,'_whole_brain_',contrast,'_0p3_',sprintf('%.02f',depth),'.mgh'];
                            if exist(projectionfile) == 0
                                projectionfile = [freesurfer_folder,'/Profiles/equi_',hemisphere,'_whole_brain_',contrast,'_0p3_',sprintf('%.01f',depth),'.mgh'];
                            end 
                            if exist(projectionfile) == 2
                                all_vals(:,d) = load_mgh(projectionfile);
                            else
                                all_vals(:,d) = NaN(length(thicknessvals),1); %%% if contrast does not exist
                            end
                        end

                        %% correct vertexwise profiles
                        for vertex = 1:size(all_vals,1)
                            profile = all_vals(vertex,:)';
                            outl = isoutlier(profile); %%% values outside the brain should be outliers
                            %%% complicated way of getting outlier indices to
                            %%% consider
                            idx = [];
                            if outl(1) == 1 %%% if it is due to bad pial segmentation, then the lowest cortical depths should be an outlier
                                idx = [idx,1];
                                allgood = 1;
                                while allgood == 1
                                   for pos = 2:length(profile)
                                       if outl(pos) == 1
                                           idx = [idx,pos];
                                       else
                                           allgood = 0;
                                       end
                                   end
                                end
                            end
                            %%% we assume that outliers are due to segmentation and
                            %%% get rid of these values, resampling the profile
                            profile_capped = profile(setdiff([1:length(profile)],idx));
                            stepsize = length(profile_capped)/1000;
                            profile_resampled = interp1(1:length(profile_capped),profile_capped,1:stepsize:length(profile_capped));
                           % stepsize = length(profile_resampled)/20;
                            stepsize = length(profile_resampled)/size(all_vals,2);
                            profile_new = interp1(1:length(profile_resampled),profile_resampled,1:stepsize:length(profile_resampled));
                            capped(vertex) = length(idx);
                            if capped(vertex) > 5 %%% something is off if too many were skipped
                                new_vals(vertex,:) = NaN * profile_new';
                            else
                                new_vals(vertex,:) = profile_new';
                            end
                        end

                        %% calculate nodewise profiles
                        for node = 1:nodenr

                            filename_profile = [outdir,'/TEST',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                            filename_profile_plot = [outdir,'/TEST',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.png'];

                            clear curr_profile curr_profile_error

                            if exist(filename_profile_plot) == 0

                             vertidx = find(label == roicodes(node));
                             if length(vertidx) > 0
                                 try
                                 subidx = find(~isoutlier(thicknessvals(vertidx)));
                                 catch
                                    1 + 1

                                 end
                                 strictidx = vertidx(subidx);

                                 %%% profile original
                                 or_pro = median(all_vals(vertidx,:))';
                                 %%% profile tidied
                                 interp_pro = nanmedian(new_vals(vertidx,:));
                                 %%% profile original thickness
                                 or_pro_thickness = median(all_vals(strictidx,:))';
                                 interp_pro_thickness = nanmedian(new_vals(strictidx,:));

                                 %%% take strictest as profile
                                 curr_profile = interp_pro_thickness;
                                 curr_profile_error = iqr(new_vals(strictidx,:));

                                   %%% exclude R1 in some
                                   if strcmp(contrast(1:2),'R1')
                                        if strcmp(brain_undersc_id(1:3), '004') || strcmp(brain_undersc_id(1:3), '032') || strcmp(brain_undersc_id(1:3), '033')
                                            curr_profile = NaN(size(curr_profile));
                                            curr_profile_error = NaN(size(curr_profile));
                                        end
                                   end

                                %%% save profile as two column file with depth and
                                %%% value
                                csvwrite(filename_profile, [depth_sampled', curr_profile', curr_profile_error']);
                                %%% plot
                                f = figure()
                                    set(f, 'visible','off');
                                    %function_plot_profile(curr_profile, curr_profile_error, depth_sampled, contrast, node)
                                    hold on
                                    plot(or_pro);
                                    plot(interp_pro);
                                    plot(or_pro_thickness);
                                    plot(interp_pro_thickness);
                                    legend('original','interpolated','original thickness','interpolated thickness');
                                    saveas(gcf,filename_profile_plot);
                                    close(f)
                             end
                         end
                    end
                 end
            end
     end
     
     %%% global profiles
     hemispheres = {'lh','rh'};
     for con = 1:length(contrasts)
            contrast = contrasts{con};
            filename_globprofile_plot = [outdir,'/TEST',brain_undersc_id,'_',atlas,'_global_',contrast,'.png'];
            filename_globprofile = [outdir,'/TEST',brain_undersc_id,'_',atlas,'_global_',contrast,'.csv'];
            contrast = contrasts{con};
            all_profiles.(contrast) = [];
            for node = 1:nodenr
                for hem = 1:2
                    hemisphere = hemispheres{hem};
                    filename_profile = [outdir,'/TEST',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
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
                    filename_profile = [outdir,'/TEST',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                    filename_partprofile = [outdir,'/TEST',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'_partial.csv'];
                    filename_partprofile_plot = [outdir,'/TEST',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'_partial.png'];
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