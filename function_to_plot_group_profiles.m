function function_to_plot_group_profiles(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, outdir_group)

   hemispheres = {'lh','rh'};
   plotcols = {'r','b','g'};
   
    for node = 0:nodenr
           %set(f, 'visible','off');
        filename_profile_plot = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_ALL.png'];   
        if exist(filename_profile_plot) == 0
        subpl_count = 1;
        for hem = 1:2
            for plottype = 1:2
               f = figure(node+1);
                set(f,'units','normalized','Position',[0 0 .9 .7]);
                subplot(2,2,subpl_count); subpl_count = subpl_count + 1;
                clear all_profiles;
                for con = 1:3%length(contrasts)
                    contrast = contrasts{con};
                    hemisphere = hemispheres{hem};
                     if plottype == 1
                        filename_group_profile = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'.csv'];
                     else
                       filename_group_profile = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_',hemisphere,'_',contrast,'_partial.csv'];
                     end
                    profdata = csvread(filename_group_profile);
                    profs_for_later(:,hem,plottype) = profdata(:,2);
                    curr_profile = profdata(:,2); curr_profile_error = profdata(:,2); depth_sampled = profdata(:,1);
                    if node > 0
                    all_profiles(:,con,node) = curr_profile;
                    end
                    hold on;
                    plot(curr_profile,depth_sampled(end:-1:1),plotcols{con},'LineWidth',2);
                    %se = shadedErrorBar(1:length(curr_profile), curr_profile, 0.5* curr_profile_error,plotcols{con},'transparent') %{'markerfacecolor','r'},0.25) 
                    %plot(curr_profile,plotcols{con},'LineWidth',2);
%                     if plottype == 1
%                         plot(curr_profile,depth_sampled(end:-1:1),plotcols{con},'LineWidth',2);
%                     %prof,rev_d
%                     %%% regress out global profile
%                     else
%                         [b bint r] = regress(curr_profile,[global_profile.(contrast),ones(length(curr_profile),1)]);
%                         plot(r,depth_sampled(end:-1:1),[plotcols{con},'--'],'LineWidth',2);
%                     end
                end
                hold off
                set(gca,'YTick',[0 max(depth_sampled)]);
                set(gca,'YTickLabel',{'WM','pial'},'FontSize',15,'FontWeight','Bold');
                xtickangle(45);
                ylabel('cortical depth','FontSize',15,'FontWeight','Bold');
                xlabel('z');
                set(gcf,'color','white');
                if node > 0
                    group_skewness = function_profile_skewness(depth_sampled', mean(all_profiles(:,:,node),2));
                    title(['node ',sprintf('%.3d',node),' sk: ', sprintf('%.04f', group_skewness)]);
                end
                legend(contrasts,'Location','southwestoutside') %'southwest')
                if plottype == 1
                   xlim([-2 2]);
                else
                   xlim([-2 2]);
                end
            end
          %  all_node_profiles(:,count) =  mean(all_profiles,2); count = count + 1;
            
        end
        %%% second plot
        1 + 1
        
        
        
        %%%
        saveas(gcf,filename_profile_plot);
        close(f)
        end
    end
    %mean(all_profiles,3)
    
end
