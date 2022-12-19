%%% to do: downsample the vertices or take a parcellation and extract
%%% values from the nodes

clear all
close all
clc

addpath(genpath('freesurfer/6.0.0/ubuntu-xenial-amd64/matlab'))
addpath(genpath('Software/BrewerMap'))

dobarplots = 0;

 brain_slash_ids = {}

hemispheres = {'lh','rh'};
contrasts = {'R1', 'MTsat','R2s'};
contrast_names = {'R1 (s^{-1})','MT_{sat} (p.u.)', 'R2* (s^{-1})'} ; 
depth_sampled = 0.05:0.05:.95;
atlas = 'BB38chimp';
nodenr = 38;
resultdir = 'Cortical_analysis';

load_chimp_ages;
load_chimp_sex;
load_chimp_wild;
load_chimp_brain_weight;

%% descriptives for entire sample
clear idx
for scan = 1:length(brain_slash_ids) % scans_to_use
    id = brain_slash_ids{scan};
    idx(scan) = str2num(id(1:3));
end
mean(ages_from_database(idx))
std(ages_from_database(idx))
sum(sex_from_database(idx) == 1)
sum(sex_from_database(idx) == 2)
sum(wild_from_database(idx) == 1)
sum(wild_from_database(idx) == 2)
sum(wild_from_database(idx) == 3)

outdir_group = [resultdir, '/group_profiles'];
example_brain_dir = [''];



%% get profiles for each brain
count = 1;
for scan = 1:length(brain_slash_ids) % scans_to_use

    brain_slash_id = brain_slash_ids{scan};
    brain_undersc_id = strrep(brain_slash_id,'/','_'); 

    id = brain_slash_ids{scan};
    idx = str2num(id(1:3));

    %%% define freesurfer dir
    freesurfer_folder = ['preprocessed/',brain_slash_id,'_V2/freesurfer/MTsat_0p3mm_downsampled_to_0p7mm_FSV7_hires']

    if exist(freesurfer_folder)
        
        %%% save profiles
        outdir_ind = [resultdir, '/individual_profiles'];
        system(['mkdir ', outdir_ind]);
        check_filename_profile = [outdir_ind,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',nodenr),'_',hemispheres{2},'_',contrasts{end},'.csv'];
        if exist(check_filename_profile) == 0
            function_to_get_depth_vals(freesurfer_folder, contrasts([1:3,5:7]), depth_sampled, atlas, nodenr, outdir_ind, brain_undersc_id)
            %function_to_get_depth_vals_test2(freesurfer_folder, contrasts([1:3]), [0:.05:1], atlas, nodenr, [outdir_ind,'_TEST'], brain_undersc_id)
        end
        
        %%% save cortical thickness
       % [thickness, surface, volume, roiname, icv, cort_vol] = function_to_get_cortical_thickness(freesurfer_folder);
        [thickness, surface, volume, normalised_volume, roiname, global_vals] = function_to_get_cortical_thickness(freesurfer_folder);
        thickness_matrix(:,count) = thickness; %%% thickness per node
        surface_matrix(:,count) = surface; %%% surface area per node
        volume_matrix(:,count) = volume; %%% volume per node
        normalised_volume_matrix(:,count) = normalised_volume; %%% in % from cortical volume
        wb.icv(count) = global_vals.ICV; %%% intracranial volume
        wb.cortical_volume(count) = global_vals.cortical_vol; %%% total cortical volume
        wb.wm_volume(count) = global_vals.wm_vol; %%% total cerebral WM volume
        
        %%%
        ages(count) = ages_from_database(idx);
        sex(count) = sex_from_database(idx);
        wild(count) = wild_from_database(idx);
        weight(count) = weight_from_database(idx);
        brain_number(count) = idx;
        count = count + 1;
    end
end

%% group thickness
meanval = mean(thickness_matrix');
stdval = std(thickness_matrix');
meanval2 = [meanval(1:38); meanval(39:end)];
stdval2 = [stdval(1:38); stdval(39:end)];
b = bar(meanval2');
hold on
er = errorbar(meanval2', stdval2');
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
er(1).Color = [0 0 0]; er(2).Color = [0 0 0];
er(1).XData = er(1).XData - 0.15; % b(1).BarWidth * 0.5;
er(2).XData = er(2).XData + 0.15; %;b(1).BarWidth * 0.5;
xticks(1:length(roiname))
xticklabels(roiname)
xtickangle(45)
hold off

%% for each ROI, contrast, get average values and average profiles (of all the ones found in the folder!!)
system(['mkdir ', outdir_group]);
depth_sampled = 0.05:0.05:0.95;
close all
%%% make profile files and plots
function_to_make_group_profiles(freesurfer_folder, contrasts([1:3,5:7]), depth_sampled, atlas, nodenr, outdir_ind, outdir_group)
%%% just for plot
function_to_make_age_grouped_profiles(freesurfer_folder, contrasts([1:3,5:7]), depth_sampled, atlas, nodenr, outdir_ind, outdir_group)
%%% load in values
[average_projection average_skewness all_table all_table_iqr all_table_skewness all_table_skewness_paquola] = function_make_average_projection(contrasts([1:3,5:7]), atlas, nodenr, outdir_group);

%%% merge tables
fn = fieldnames(all_table_skewness);
for v = 1:length(fn)
   newname = ['skewness_',fn{v}];
   all_table.(newname) = all_table_skewness.(fn{v});
   newname = ['skewness_paquola_',fn{v}];
   all_table.(newname) = all_table_skewness_paquola.(fn{v});
end
% 
% %%% quick test
% for i = 1:76
%     scatter(all_table.skewness_MTsat(:,i), all_table.skewness_paquola_MTsat(:,i));
%     lsline
%     pause(1)
% end

globproffile = [outdir_group,'/globalprofiles.mat'];
if exist(globproffile) == 2
    load(globproffile)
else
    globalprofile = function_make_global_profile(contrasts([1:3,5:7]), atlas, nodenr, outdir_group)
    save(globproffile,'globalprofile')
end

%%% save global as node 00
types = {'','_partial'}
for t = 1:2
     type = types{t};
    for c = [1:3, 5:7]         
        contrast = contrasts{c}
        for hem = 1:2
            hemisphere = hemispheres{hem};
            filename_group_profile = [outdir_group,'/Group_profiles_',atlas,'_000_',hemisphere,'_',contrast,type,'.csv'];
           %%% write with error0 
            csvwrite(filename_group_profile, [depth_sampled',globalprofile.(contrast), zeros(length(globalprofile.(contrast)),1)]);
        end
    end
end
%%%
%function_to_plot_group_profiles(freesurfer_folder, contrasts, depth_sampled, atlas, nodenr, outdir_group)
function_to_plot_group_profiles_average_LR(freesurfer_folder, contrasts([1:3,5:7]), depth_sampled, atlas, nodenr, outdir_group)
function_to_plot_group_profiles_average_extended(freesurfer_folder, contrasts([1:3,5:7]), depth_sampled, atlas, nodenr, outdir_group)

%%% make matrix for plotting 
contrasts{9} = 'thickness'; 
contrast_names{9} = 'thickness (mm)';
all_table.thickness = transpose(thickness_matrix);
contrasts{10} = 'surface'
contrast_names{10} = 'surface (mm2)';
all_table.surface = transpose(surface_matrix);
contrasts{11} = 'volume'
contrast_names{11} = 'volume (mm3)';
all_table.volume = transpose(volume_matrix);
contrasts{12} = 'normalised volume'
contrast_names{12} = 'normalised volume (mm3)';
all_table.normalised_volume = transpose(normalised_volume_matrix);


%% deal with whole brain measures
subplot(1,3,1)
    scatter(weight, wb.icv)
    xlabel('weight in g')
    ylabel('intracranial volume from freesurfer')
    lsline
    hold on
    text(weight, wb.icv, cellstr(num2str(brain_number')))
subplot(1,3,2)
    scatter(weight, wb.cortical_volume)
    xlabel('weight in g')
    ylabel('cortical volume from freesurfer')
    lsline
    hold on
    text(weight, wb.cortical_volume, cellstr(num2str(brain_number')))  
subplot(1,3,3)
    scatter(wb.wm_volume, wb.cortical_volume)
    xlabel('wm volume from freesurfer')
    ylabel('cortical volume from freesurfer')
    lsline
    hold on
    text(wb.wm_volume, wb.cortical_volume, cellstr(num2str(brain_number')))  
   saveas(gcf,[resultdir,'/figures/Wholebrain_measures_scatters.png']);    

%% model and plot
volume_measures = fieldnames(wb);
figure()
for v = 1:length(volume_measures)
    subplot(1,length(volume_measures),v);
    volume_measure = volume_measures{v};
    iv = ages;
    dv = wb.(volume_measure);
    %%% this part of the code is copied from below
    fake_iv = 0:0.01:max(iv);
    %%% hallgren model as it is
    newfun = @(cpars,iv)(cpars(1)*(1-exp(cpars(2)*iv*-1))+cpars(3));
    pars0 = [5, 0.1, 1]; %
    [BETA,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(newfun,pars0,iv,dv);
    ci = nlparci(BETA, RESIDUAL,'Jacobian',JACOBIAN);
    prediction = newfun(BETA, iv);
    fake_pred = newfun(BETA, fake_iv);
    scatter(iv,dv,20,'k','filled')
    hold on
    aah = plot(fake_iv,fake_pred,'Color','k','LineWidth',3.5);
   %lsline
   ylabel(volume_measure,'interpreter','none');
    xlabel('age(years)');
    %aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
 saveas(gcf,[resultdir,'/figures/Lineplot_wholebrain_measures.png']);  



%% now ROI based analysis

mask_idx = 1:38;

%%% do babies
babyfile = [resultdir,'/BabyROIvals_median.mat'];
if exist(babyfile) == 0
    count2 = 1;
    for scan2 = [2 9]

        brain_slash_id = brain_slash_ids{scan2};
        brain_undersc_id = strrep(brain_slash_id,'/','_'); 

        id = brain_slash_ids{scan2};
        idx = str2num(id(1:3));
        ages2(count2) = ages_from_database(idx);

        ROI_folder = ['preprocessed/',brain_slash_id,'_V2/manual_ROIs'];
        MPM_folder = ['preprocessed/',brain_slash_id,'_V2/MPMs_to_use'];
        for c = 1:4
            contrast = contrasts{c};
            %%% constant value across brain (robust, but ignores variation)
             %%% add babies
            file = [MPM_folder,'/',contrast,'_0p3_run01_brain_masked_reoriented.nii']
            %age_vector_appended = [age_vector; ages_from_database(2)];
            [error val] = system(['fslstats ', file, ' -P 50']);
            [error valu] = system(['fslstats ', file, ' -P 75']);
            [error vall] = system(['fslstats ', file, ' -P 25']);
            all_baby_table.(contrast)(:,count2) = repmat(str2double(val),nodenr*2,1);
            all_baby_iqr_table.(contrast)(:,count2) = repmat(str2double(valu)-str2double(vall),nodenr*2,1);
        end
        count2 = count2 + 1;
    end
    %all_baby_table = function_to_get_
    save(babyfile,'all_baby_table','all_baby_iqr_table','ages2');
else
    load(babyfile,'all_baby_table','all_baby_iqr_table','ages2');
end

%% put all values in tables
all_variables = fieldnames(all_table);

cs = 1:length(all_variables);
value_name = 'test';
clear vals vals_left vals_right vals_mean vals_iqr_left vals_iqr_right vals_iqr_mean
for c = cs
    contrast = all_variables{c};
       for region = 1:length(mask_idx)
            %curr_con = contrasts{c};
            origvals = all_table.(contrast); 
            %%% median
            try
                vals_baby = all_baby_table.(contrast)';
                origvals = [origvals; vals_baby];
                [ages_sorted, sort_idx] = sort([ages, ages2]);
            catch
                [ages_sorted, sort_idx] = sort([ages]); 
            end
            vals.left.(contrast)(:,region) = origvals(sort_idx,mask_idx(region));
            vals.right.(contrast)(:,region) = origvals(sort_idx,mask_idx(region)+38);
            %%% Iqr
            try
                vals_iqr = all_table_iqr.(contrast); 
                vals_baby_iqr = all_baby_iqr_table.(contrast)';
                vals_iqr = [vals_iqr; vals_baby_iqr];
                [ages_sorted, sort_idx] = sort([ages, ages2]);
            catch
                vals_iqr = zeros(size(origvals));
                [ages_sorted, sort_idx] = sort([ages]);
            end
            vals.iqr_left.(contrast)(:,region) = vals_iqr(sort_idx,mask_idx(region));
            vals.iqr_right.(contrast)(:,region) = vals_iqr(sort_idx,mask_idx(region)+38);
            vals.age.(contrast) = ages_sorted';
       end
       vals.mean.(contrast) = nanmean(cat(3 ,vals.left.(contrast), vals.right.(contrast)), 3); 
       vals.iqr_mean.(contrast) = nanmean(cat(3 ,vals.iqr_left.(contrast), vals.iqr_right.(contrast)), 3); 
end        

%% fit;
clear model model2
hemispheres = {'left','right','mean'}
for region = 1:38 %length(mask_idx)
    for h = 1:length(hemispheres)
        hemi = hemispheres{h};
         for c = cs 
            contrast = all_variables{c};
            values_to_use = vals.(hemi).(contrast)(:,region);
            
            
            %% test a model
            iv = vals.age.(contrast);
            dv = values_to_use;
            
            fake_iv = 0:0.01:max(iv);
            
            %%% hallgren model as it is
            newfun = @(cpars,iv)(cpars(1)*(1-exp(cpars(2)*iv*-1))+cpars(3));
            pars0 = [5, 0.1, 1]; %%% parameter initialisation based on rough guess of paper, 1 = amplitude, 3 = constant, 2 = rate

            if strcmp(contrast,'thickness') 
                 %%% adapted model for thickness
                 newfun = @(cpars,iv)(cpars(1)*(exp(cpars(2)*iv*-1))+cpars(3)); %%% exponential decay
            end
            warning('') %% stop if there is warning
            %% fit linear model
            if strcmp(contrast,'surface') || strcmp(contrast,'volume') || strcmp(contrast(1:2),'sk')   %%% model with a linear model 
                %pars0 = [0 0 0];
                %newfun = @(cpars,iv)(0*cpars(3) + cpars(2)*iv); %%% just linear model
                %%% i want rate to be the second and constant the
                %%% 3rd parameter, but it fails with 3 parameter
                %%% initialisations, so run with 2 and add 0 to
                %%% vecor later
                pars0 = [0 0];
                newfun = @(cpars,iv)(cpars(2) + cpars(1)*iv); %%% i want rate be second and constant 3rd parameter
                %[BETA,R,J,COVB,MSE] = nlinfit(iv,dv,newfun,pars0); 
                %%% use different fitting function that can constrain
                %%% parameters; do constrain constant to be positive for surface and volume
                nonnanidx = find(~isnan(dv));
                if strcmp(contrast,'surface') || strcmp(contrast,'volume') %%% lower bound for constant = 0
                   lb = [-9999, 0]; ub = [9999 9999];
                   [BETA,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(newfun,pars0,iv(nonnanidx),dv(nonnanidx),lb,ub);
                elseif strcmp(contrast(1:2),'sk') %%% skewness no constraints
                   %lb = [-9999, -9999]; ub = [9999 0]; %%% constrain to negative
                   lb = [-9999, -9999]; ub = [9999 9999]; %%% no constrain
                   [BETA,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(newfun,pars0,iv(nonnanidx),dv(nonnanidx));
                end 
                %BETA
                %X
                %ci = nlparci(BETA, R,'Jacobian',J) %%% confidence intervals of parmaeters
                ci = nlparci(BETA, RESIDUAL,'Jacobian',JACOBIAN);
                prediction = newfun(BETA, iv);
                fakeprediction = newfun(BETA, fake_iv);
                %%% add fake amplitude 
                BETA = [0 BETA]; %%% make rate second and constant thrid
                ci = [0 0; ci];
            else
                %% exponential model
                try
                    %[BETA,R,J,COVB,MSE] = nlinfit(iv,dv,newfun,pars0); 
                    %ci = nlparci(BETA, R,'Jacobian',J); %%% confidence intervals of parmaeters
                    nonnanidx = find(~isnan(dv));
                    
                    lb = [-9999, -9999, -9999]; ub = [9999, 9999,  9999]; %%% constrain constant to be positive does not help because the amplitude parameter is also important
                    %%% constrain rate (second parameter) to be positive
                    if strcmp(contrast(1:2),'th') %%% thickness cannot be more than 5 mm
                        lb = [-9999, 0, 0]; ub = [5, 9999, 5]; %%% constrain constant to be positive does not help because the amplitude parameter is also important
                    elseif strcmp(contrast(1:2),'sk') %%% skewness does not need constraitnss
                        lb = [-10, -9999, -10]; ub = [10, 9999, 10]; %%% skewness is theoretically negative
                    else %%% for MPMs constrain constant and amplitude to min of 0
                        lb = [0, 0, 0]; ub = [9999, 9999, 9999];
                    end
                    
                    [BETA,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(newfun,pars0,iv(nonnanidx),dv(nonnanidx),lb,ub);
                    ci = nlparci(BETA, RESIDUAL,'Jacobian',JACOBIAN);
                    prediction = newfun(BETA, iv);
                    fakeprediction = newfun(BETA, fake_iv);
                catch
                    BETA = NaN(1,length(pars0));
                    ci = NaN(length(pars0),2);
                    prediction = NaN(size(iv));
                    fakeprediction = NaN(size(fake_iv));
                end
            end
           
           model.(contrast).(hemi).rate_CI_low(region) = ci(2,1); %%% steepness of incline
            model.(contrast).(hemi) = model.(contrast).(hemi);
           model.(contrast).(hemi).rate_CI_high(region) = ci(2,2);    
            model.(contrast).(hemi).amplitude_CI_low(region) = ci(1,1); %%% amplitude
            model.(contrast).(hemi).amplitude_CI_high(region) = ci(1,2);       
            model.(contrast).(hemi).constant_CI_low(region) = ci(3,1); %%% arbitrary constant
            model.(contrast).(hemi).constant_CI_high(region) = ci(3,2); 

            
            [r p] = corrcoef(dv,prediction,'rows','complete');
            model.(contrast).(hemi).coeffr2(region) = r(1,2).^2;
            model.(contrast).(hemi).p(region) = p(1,2);
            model.(contrast).(hemi).n(region) = sum(~isnan(dv));
            
            model.(contrast).(hemi).rate(region) = BETA(2);
            model.(contrast).(hemi).amplitude(region) = BETA(1);
            model.(contrast).(hemi).constant(region) = BETA(3);

             [warnMsg, warnId] = lastwarn;
             if ~isempty(warnMsg)
                1 + 1 
                %[contrast,' ',value_name,' ',mask_names2{region}]
                model.(contrast).(hemi).problem(region) = 1;
             else
                model.(contrast).(hemi).problem(region) = 0;
             end
            
             model.(contrast).(hemi).iv(:,region) = iv;
             model.(contrast).(hemi).dv(:,region) = dv;
             model.(contrast).(hemi).pred(:,region) = prediction;
             model.(contrast).(hemi).fake_pred(:,region) = fakeprediction;
             model.(contrast).(hemi).fake_iv(:,region) = fake_iv;
             
             %% compare to Miller model
            %%% iv = age
            iv_trunc = iv;
            iv_trunc(iv_trunc > 17) = 17;
            fake_iv = 0:0.01:max(iv_trunc);
            ivmat = [ones(length(iv_trunc),1), iv_trunc, iv_trunc.^2];
            [B,BINT,R,RINT,STATS] = regress(dv,ivmat);
            predict = ivmat * B;
            fakepred = [ones(length(fake_iv),1), fake_iv', fake_iv'.^2] * B;
            
            model2.(contrast).(hemi).iv(:,region) = iv_trunc;
            model2.(contrast).(hemi) = model2.(contrast).(hemi);
            model2.(contrast).(hemi).dv(:,region) = dv;
            model2.(contrast).(hemi).pred(:,region) = predict;
            model2.(contrast).(hemi).fakepred(:,region) = fakepred;
            model2.(contrast).(hemi).fake_iv(:,region) = fake_iv;
            model2.(contrast).(hemi).r2(:,region) = STATS(1);
            model2.(contrast).(hemi).p(:,region) = STATS(3);

  
         end
    end
end

%% plot for selected regions
mask_names2 = {'Visual','Auditory','Motor','Somatosensory','Broca homologue','Prefrontal'};
mask_idx2 = [24 34 1 26 5 8]; 
six_color_codes = {[140 20 60], [26 58 148], [60 140 180], [20 180 140], [220 20 100], [180 22 140]};


%% line plots
hemi = 'mean'
for c = [1:3, 7, 9, 11, 19] %1:3 %cs %s = 1:length(all_variables)
    f1 = figure();
     set(f1,'Position',[100 100 1900 700]);
    f2 = figure(); 
    set(f2,'Position',[100 100 1500 700]);
    for regiontoplot = 1:length(mask_idx2)
        contrast = all_variables{c};
        region_idx = mask_idx2(regiontoplot);
        
        thismodel = model.(contrast).(hemi);
        iv = thismodel.iv(:,region_idx);
        dv = thismodel.dv(:,region_idx);
        fake_iv = thismodel.fake_iv(:,region_idx);
        fake_pred = thismodel.fake_pred(:,region_idx);  
        pred = thismodel.pred(:,region_idx);
        rate = thismodel.rate(:,region_idx);
        constant = thismodel.constant(:,region_idx);
        amplitude = thismodel.amplitude(:,region_idx);
        pval = thismodel.p(:,region_idx);
        r2 = thismodel.coeffr2(:,region_idx);
        
        %% plot hallgren model
        figure(f1)
        hold on
        scatter(iv, dv, 100, six_color_codes{regiontoplot}./255,'filled');
        if pval < 0.05
          if rate ~= 0
              if strcmp(contrast,'thickness') 
                 legend_entry{regiontoplot} = [mask_names2{regiontoplot},': constant: ',num2str(round(constant,2)),', time constant: ' ,num2str(round(1/rate,2))];
              elseif strcmp(contrast,'surface') || strcmp(contrast,'volume') || strcmp(contrast(1:2),'sk')
                legend_entry{regiontoplot} = [mask_names2{regiontoplot},': ', num2str(round(constant,2)),' + age x ' , num2str(round(rate,4))];
              else
                 legend_entry{regiontoplot} = [mask_names2{regiontoplot},': amplitude: ',num2str(round(amplitude + constant,2)),', time constant: ' ,num2str(round(1/rate,2))];
              end 
         else
             legend_entry{regiontoplot} = [mask_names2{regiontoplot}];
          end
        else
           legend_entry{regiontoplot} = [mask_names2{regiontoplot},': n.s.'];
        end
        aah = plot(fake_iv,fake_pred,'Color',six_color_codes{regiontoplot}./255,'LineWidth',3.5);
        aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold off
        xlabel('years');
        if strcmp(contrast,'R1');
            conlab = 'R1 (s^{-1})';
        elseif strcmp(contrast,'R2s');
            conlab = 'R2* (s^{-1})';
        elseif strcmp(contrast,'MTsat');
            conlab = 'MT_{sat} (p.u.)';
            conlab = 'MTsat (p.u.)';
        elseif strcmp(contrast,'skewness_R1');
            conlab = 'profile skewness R1';
        elseif strcmp(contrast,'skewness_R2s');
            conlab = 'profile skewness R2*';   
        elseif strcmp(contrast,'skewness_MTsat');
            conlab = 'profile skewness MT_{sat}';
            conlab = 'profile skewness MTsat'; 
        else
            conlab = contrast;
        end
        ylabel(conlab);
        set(gca, 'FontSize',35);
        set(gca, 'Color', 'white');
        all_min(regiontoplot) = min([dv]);
        all_max(regiontoplot) = max([dv]);
        if min(all_min) > 0 %%% only if no negative values
            try
                m = max(all_max);
                ylim([0, max(m) + 0.2*max(m)]);           
            end
        end
        
        %% miller model
        thismodel2 = model2.(contrast).(hemi);
        iv = thismodel2.iv(:,region_idx);
        dv = thismodel2.dv(:,region_idx);
        pred = thismodel2.pred(:,region_idx);
        fakepred = thismodel2.fakepred(:,region_idx);
        fake_iv = thismodel2.fake_iv(:,region_idx);
        r2val = thismodel2.r2(:,region_idx);
        figure(f2)
            hold on;
            scatter(iv,dv, 100, six_color_codes{regiontoplot}./255,'filled');
            h = plot(fake_iv,fakepred,'Color',six_color_codes{regiontoplot}./255,'LineWidth',3.5);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            hold off
            xlabel('years');
            ylabel(contrast, 'Interpreter', 'none');
            set(gca, 'FontSize',25);
            set(gca, 'Color', 'white');
            legend_entry2{regiontoplot} = [mask_names2{regiontoplot},': ', num2str(round(r2val,2))]; 
        
    end
 
    figure(f1)
        legend(legend_entry,'Location','NorthWestOutside','FontSize',15);
        saveas(gcf,[resultdir,'/figures/Lineplot_',contrast,'_',hemi,'.png']);
    figure(f2)
        legend(legend_entry2,'Location','NorthWestOutside','FontSize',15);
        saveas(gcf,[resultdir,'/figures/Lineplot_QuadraticFit_',contrast,'_',hemi,'.png']);  

end

%% for exploration, make trajectory plots for individual regions
hemi = 'mean'
for c = 1:3 
    contrast = all_variables{c};
    f1 = figure();
     set(f1,'Position',[100 100 1900 700]);
     subplcount = 1;
    for region_idx = 1:38
        thismodel = model.(contrast).(hemi);
        iv = thismodel.iv(:,region_idx);
        dv = thismodel.dv(:,region_idx);
        fake_iv = thismodel.fake_iv(:,region_idx);
        fake_pred = thismodel.fake_pred(:,region_idx);  
        pred = thismodel.pred(:,region_idx);
        rate = thismodel.rate(:,region_idx);
        constant = thismodel.constant(:,region_idx);
        amplitude = thismodel.amplitude(:,region_idx);
        pval = thismodel.p(:,region_idx);
        r2 = thismodel.coeffr2(:,region_idx);
        
        %% plot hallgren model
        figure(f1)
        subplot(6,7,subplcount); subplcount = subplcount + 1;
        hold on
        scatter(iv, dv, 30, six_color_codes{regiontoplot}./255,'filled');
        aah = plot(fake_iv,fake_pred,'Color',six_color_codes{regiontoplot}./255,'LineWidth',1.5);
        aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold off
        xlabel('years');
        if strcmp(contrast,'R1');
            conlab = 'R1 (s^{-1})';
        elseif strcmp(contrast,'R2s');
            conlab = 'R2* (s^{-1})';
        elseif strcmp(contrast,'MTsat');
            conlab = 'MTsat (p.u.)';
        else
            conlab = contrast;
        end
        ylabel(conlab);
        set(gca, 'FontSize',10);
        set(gca, 'Color', 'white');
        all_min(regiontoplot) = min([dv]);
        all_max(regiontoplot) = max([dv]);
        if min(all_min) > 0 %%% only if no negative values
            try
                m = max(all_max);
                ylim([0, max(m) + 0.2*max(m)]);           
            end
        end        
    end
end

%% make output table
%%% to check p vals for all regions:
table_types = {'6rois','all'};
for tabtodo = 1:2
    if strcmp(table_types{tabtodo}, 'all')
        mask_idx2 = 1:38;
        mask_names2 = num2cell([1:38]);
        annotfile = ['results/Cortical_analysis/group_freesurfer/fsaverage/label/lh.BB38chimp.annot']
        [vertices, labeling, colortable] = read_annotation(annotfile);
        mask_names2 = colortable.struct_names(2:end)';
    else
        %%% for paper table
        mask_names2 = {'Visual','Auditory','Motor','Somatosensory','Broca homologue','Prefrontal'};
        mask_idx2 = [24 34 1 26 5 8]; 
    end
    firstcol = mask_names2';
    mat1 = [];
    for c = [1:3, 7, 9, 11, 19] %cs
     contrast = all_variables{c}
     shortcut = model.(contrast).(hemi);
     miller_r2 = model2.(contrast).(hemi).r2;
     perc_at_birth = 100 * shortcut.constant(mask_idx2) ./ (shortcut.constant(mask_idx2) + shortcut.amplitude(mask_idx2));
     mat1 = [shortcut.n(mask_idx2)',round(shortcut.coeffr2(mask_idx2),2)',round(shortcut.p(mask_idx2),3)',round(shortcut.constant(mask_idx2),2)',round(shortcut.amplitude(mask_idx2),2)', round(shortcut.constant(mask_idx2) + shortcut.amplitude(mask_idx2),2)', round(perc_at_birth)',round(1./shortcut.rate(mask_idx2),2)',round(shortcut.rate(mask_idx2),4)',round(shortcut.rate_CI_low(mask_idx2),2)',round(shortcut.rate_CI_high(mask_idx2),2)',round(miller_r2(mask_idx2)',2)];
     table_to_write = cell2table([firstcol,num2cell(mat1)],'VariableNames',{'Region','n','R^2', 'p', 'constant', 'amplitude', 'plateau', 'perc', 'time const', 'rate', 'lower rate CI', 'higher rate CI','MillerR2'})
     writetable(table_to_write,[resultdir,'/model_results_',contrast,'_',hemi,'_',table_types{tabtodo}]);
    end  
end



%%% average time constant values across cortex for mtsag
median(1./model.MTsat.mean.rate)
3.1028 * 3 %%% highest realistic value x 3


%%% relationship between rate and amplitude
varstotake = {'R2s','MTsat'};
for v = 1:2
    f = figure()
    set(f,'Position',[100 100 300 300])
    var = varstotake{v};
    pvals = [model.(var).left.p, model.(var).right.p];
    a1 = [1./model.(var).left.rate, 1./model.(var).right.rate]; %%% rate
    a2 = [model.(var).left.rate, model.(var).right.rate]; %%% constant
    b = [model.(var).left.constant + model.(var).left.amplitude, model.(var).right.constant + model.(var).right.amplitude];
    %%%
    idx = find(pvals < .05 & ~isoutlier(a1) & ~isoutlier(b)); %%% remove unrealistic values
    %idx = find(pvals < .05);
    figure()
    scatter(a2(idx),b(idx),'k','filled');
    [rv pv] = corrcoef(a2(idx),b(idx));
    xlabel('Rate');
    ylabel('amplitude');
    title([var, ' r = ', num2str(round(rv(1,2),2)), '; p = ', num2str(round(pv(1,2),3))]);
    lsline
      set(gca, 'FontSize',25);
    set(gca, 'Color', 'white'); 
    saveas(gcf,[resultdir,'/figures/Scatter_',var,'_rate_vs_plat.png']);
    %%
    figure()
    idx = find(pvals < .05 & ~isoutlier(a2) & ~isoutlier(b));
    scatter(a1(idx),b(idx),'k','filled');
    [rv pv] = corrcoef(a1(idx),b(idx));
    xlabel('Time constant'); 
    ylabel('amplitude');
    title([var, ' r = ', num2str(round(rv(1,2),2)), '; p = ', num2str(round(pv(1,2),3))]);
    lsline
    set(gca, 'FontSize',25);
    set(gca, 'Color', 'white');
    saveas(gcf,[resultdir,'/figures/Scatter_',var,'_timeconst_vs_plat.png']);
end


%% check meaningfulness of parameters by correlating left with right
variables = {'MTsat', 'R2s', 'surface', 'volume', 'thickness', 'skewness_R2s'};
for v = 1:length(variables)
    subplot(3,2,v)
    variable = variables{v};
    scatter(model.(variable).left.rate, model.(variable).right.rate);
    lsline
    xlabel('rate left');
    ylabel('rate right');
    title(variable)
end
saveas(gcf,[resultdir,'/figures/Scatter_model_pars_left_vs_right.png']);


variablestruc = {[1./model.R2s.left.rate, 1./model.R2s.right.rate], [model.R2s.left.rate, model.R2s.right.rate],[model.R2s.left.constant + model.R2s.left.amplitude, model.R2s.right.constant + model.R2s.right.amplitude] ...
    [1./model.MTsat.left.rate, 1./model.MTsat.right.rate],[model.MTsat.left.rate, model.MTsat.right.rate],  [model.MTsat.left.constant + model.MTsat.left.amplitude, model.MTsat.right.constant + model.MTsat.right.amplitude] ...
    median(all_table_skewness.R2s), median(all_table_skewness.MTsat),nanmedian(all_table_skewness.R1),  ontogenetic_surface_expansion_perc}; %, phylogenetic_surface_expansion};
variable_names = {'R2s_time_constant','R2s_rate','R2s_amplitude','MTsat_time_constant','MTsat_rate','MTsat_amplitude','R2s_skewness','MTsat_skewness','R1_skewness','ontogenetic_surface_expansion'}; %,'phylogenetic_surface_expansion' };



%%% run in terminal in FS environment to get annotatoins on template
%['mris_ca_label -t ',atlas_folder,'/BB38chimp.annot.ctab fsaverage lh ',freesurfer_folder,'/surf/lh.sphere.reg ', atlas_folder,'/lh.BB38chimp.gcs ',freesurfer_folder,'/label/lh.BB38chimp.annot']
%['mris_ca_label -t ',atlas_folder,'/BB38chimp.annot.ctab fsaverage rh ',freesurfer_folder,'/surf/rh.sphere.reg ', atlas_folder,'/rh.BB38chimp.gcs ',freesurfer_folder,'/label/rh.BB38chimp.annot']

for var = 1:length(variable_names)
    
    variable_name = variable_names{var};
    variable_vec = variablestruc{var};
    figure_file_to_save_prefix = [resultdir,'/roi_wise_projection_figures/fsaverage_',hemisphere,'_',variable_name,'_brewermap'];
    projection_file_to_save = [resultdir,'/roi_wise_projection_figures/fsaverage_',hemisphere,'_',variable_name,'.mgh'];
    
    function_make_roiwise_surface_plot(variable_vec, variable_name, [1 99], projection_file_to_save, figure_file_to_save_prefix)

end



close all
save([resultdir,'/microstructure_pipeline_workspace'])
%adult_surf_surface_matrix(:,ages>17)
