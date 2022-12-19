function function_to_plot_group_profiles_average_LR(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, outdir_group)

   hemispheres = {'lh','rh'};
   plotcols = {'r','b','g'};
   plotcols = {[123 15 15]./255, [20 15 123]./255, [19 92 45]./255};
   
    for node = 0:nodenr
        %set(f, 'visible','off');
        filename_profile_plot = [outdir_group,'/Group_profiles_',atlas,'_',sprintf('%.3d',node),'_ALL.png'];   
        if exist(filename_profile_plot) == 0
        %%% get data
        clear profs_for_later profs_error_for_later profs_depths_sampled_for_later
        for hem = 1:2
            for plottype = 1:2
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
                    profs_for_later(:,con,hem,plottype) = curr_profile;
                    profs_error_for_later(:,con,hem,plottype) = curr_profile_error;
                    profs_depths_sampled_for_later(:,con,hem,plottype) = depth_sampled;
                end
            end
        end
        f = figure(node+1);
        set(f,'Position',[100 100 1200 300]);
        %set(f,'units','normalized','Position',[0 0 .2 .45]);
        subpl_count = 1;
        for plottype = 1:2  
            subplot(1,2,subpl_count); subpl_count = subpl_count + 1;
            for con = 1:3
                 curr_profile = mean(squeeze(profs_for_later(:,con,:,plottype)),2); %%% average LH and RH
                 curr_profile_error = mean(squeeze(profs_error_for_later(:,con,:,plottype)),2); 
                 depth_sampled = profs_depths_sampled_for_later(:,con,1,plottype);
                 
                 hold on;
                 plot(curr_profile,depth_sampled(end:-1:1),'Color',plotcols{con},'LineWidth',2);
                 %se = shadedErrorBar(1:length(curr_profile), curr_profile, 0.5* curr_profile_error,plotcols{con},'transparent') 
                 %se = shadedErrorBar(curr_profile, depth_sampled(end:-1:1), 0.5* curr_profile_error,plotcols{con},'transparent') 
                 %plot(curr_profile,depth_sampled(end:-1:1),plotcols{con},'LineWidth',2);
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
            set(gca,'YTick',[0.05 0.25 0.5 0.75 max(depth_sampled)]);
            set(gca,'YTickLabel',{'WM','25%','50%','75%','pial'},'FontSize',16); %,'FontWeight','Bold');
            xtickangle(45);
            ylabel('cortical depth','FontSize',16) ; %,'FontWeight','Bold');
            xlabel('z');
            set(gcf,'color','white');
            legend(contrasts(2:3),'Location','southwestoutside') %'southwest')
            xlim([-2 2]);           
        end
        saveas(gcf,filename_profile_plot);
        close(f)
        end
    end
end
