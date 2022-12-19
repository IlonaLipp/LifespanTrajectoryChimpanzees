function function_plot_profile(curr_profile, curr_profile_error, depth_sampled, contrast, node)
    hold on;
    se = shadedErrorBar(1:length(curr_profile), curr_profile, 0.5* curr_profile_error,'b','transparent') %{'markerfacecolor','r'},0.25) 
    plot(curr_profile,'b-');
    hold off;
    set(gca,'XTick',[1,length(depth_sampled)]);
    set(gca,'XTickLabel',{'pial','WM'},'FontSize',15,'FontWeight','Bold');
    xtickangle(45);
    xlabel('cortical depth','FontSize',15,'FontWeight','Bold');
    ylabel(contrast);
    title(['node ',sprintf('%.3d',node)]);
    set(gcf,'color','white');
end