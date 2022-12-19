function global_profile = function_make_global_profile(contrasts, atlas, nodenr, outdir_group)

%%% averages all profiles from all nodes and hemispheres

   hemispheres = {'lh','rh'};
   types = {'','_partial'}
   for t = 1
     type = types{t};
     for con = 1:length(contrasts)
            contrast = contrasts{con};
            profcount = 0;
            count = 1;
            for node = 1:nodenr
                
                for hem = 1:2
                      
                    hemisphere = hemispheres{hem};
                    filename_group_profile = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,type,'.csv'];
                    datamat = importdata(filename_group_profile);
                    all_profiles(:,count) = datamat(:,2); count = count + 1;
                end
            end
            global_profile.(contrast) = mean(all_profiles,2);
     end
   end
   
end