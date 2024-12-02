%% Microstructural Application Pipeline
% This script processes and analyzes microstructural brain imaging data, particularly focusing on
% developmental trajectories across different brain regions.
%
% Key Features:
% - Processes multiple brain imaging contrasts (R1, MTsat, R2s)
% - Analyzes cortical thickness, surface area, and volume measurements
% - Fits developmental trajectories using exponential and linear models
% - Compares left and right hemisphere measurements
% - Generates various visualizations and statistical analyses
%
% Dependencies:
% - FreeSurfer (v6.0.0)
% - BrewerMap
% - Custom functions for processing and analysis

% Initialize environment
clear all
close all
clc

% Add required paths
addpath(fullfile(pwd, '../functions'));
addpath(genpath('freesurfer/6.0.0/ubuntu-xenial-amd64/matlab'))
addpath(genpath('Software/BrewerMap'))

% Configuration
dobarplots = 0;
brain_slash_ids = {} % List of brain IDs to process
hemispheres = {'lh','rh'};
contrasts = {'R1', 'MTsat','R2s'};
contrast_names = {'R1 (s^{-1})','MT_{sat} (p.u.)', 'R2* (s^{-1})'};
depth_sampled = 0.05:0.05:.95;
atlas = 'BB38chimp';
nodenr = 38;
resultdir = 'Cortical_analysis';

% Load demographic and metadata from files
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process brain scans
% Iterate through brains to:
% - Generate depth profiles if needed
% - Extract cortical measurements
% - Store regional and global metrics
% - Record demographics


%% get profiles for each brain
% Initialize counter for valid brain scans
count = 1;

% Iterate through all brain scan IDs
for scan = 1:length(brain_slash_ids) % scans_to_use

    % Get brain ID and convert slashes to underscores for filenames
    brain_slash_id = brain_slash_ids{scan};
    brain_undersc_id = strrep(brain_slash_id,'/','_'); 

    % Extract numeric ID from brain scan name
    id = brain_slash_ids{scan};
    idx = str2num(id(1:3));

    % Define FreeSurfer directory path for this brain
    freesurfer_folder = ['preprocessed/',brain_slash_id,'_V2/freesurfer/MTsat_0p3mm_downsampled_to_0p7mm_FSV7_hires']

    % Only process if FreeSurfer folder exists
    if exist(freesurfer_folder)
        
        % Create directory for individual profiles if needed
        outdir_ind = [resultdir, '/individual_profiles'];
        system(['mkdir ', outdir_ind]);

        % Check if profile already exists, generate if not
        check_filename_profile = [outdir_ind,'/',brain_undersc_id,'_',atlas,'_',sprintf('%.3d',nodenr),'_',hemispheres{2},'_',contrasts{end},'.csv'];
        if exist(check_filename_profile) == 0
            function_to_get_depth_vals(freesurfer_folder, contrasts([1:3,5:7]), depth_sampled, atlas, nodenr, outdir_ind, brain_undersc_id)
        end
        
        % Extract cortical measurements from FreeSurfer data
        [thickness, surface, volume, normalised_volume, roiname, global_vals] = function_to_get_cortical_thickness(freesurfer_folder);
        
        % Store regional measurements in matrices
        thickness_matrix(:,count) = thickness;      % Thickness per node
        surface_matrix(:,count) = surface;          % Surface area per node 
        volume_matrix(:,count) = volume;            % Volume per node
        normalised_volume_matrix(:,count) = normalised_volume; % Volume as % of cortical volume
        
        % Store global brain measurements
        wb.icv(count) = global_vals.ICV;                    % Intracranial volume
        wb.cortical_volume(count) = global_vals.cortical_vol; % Total cortical volume
        wb.wm_volume(count) = global_vals.wm_vol;           % Total cerebral white matter volume
        
        % Store demographic and metadata
        ages(count) = ages_from_database(idx);      % Age
        sex(count) = sex_from_database(idx);        % Sex
        wild(count) = wild_from_database(idx);      % Wild status
        weight(count) = weight_from_database(idx);  % Brain weight
        brain_number(count) = idx;                  % Brain ID number
        
        % Increment counter for next valid brain
        count = count + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot mean cortical thickness by region
% Create bar plot showing mean cortical thickness across brain regions,
% with error bars representing standard deviation. Data is split into
% hemispheres and formatted with angled ROI labels for readability.

%% group thickness
% Calculate mean and standard deviation of thickness across all brains
meanval = mean(thickness_matrix');
stdval = std(thickness_matrix');

% Split means and std devs into two groups (first 38 and remaining ROIs)
meanval2 = [meanval(1:38); meanval(39:end)];
stdval2 = [stdval(1:38); stdval(39:end)];

% Create bar plot of mean thicknesses
b = bar(meanval2');
hold on

% Add error bars showing standard deviation
er = errorbar(meanval2', stdval2');

% Format error bars - make them black with no connecting lines
er(1).LineStyle = 'none'; er(2).LineStyle = 'none'; 
er(1).Color = [0 0 0]; er(2).Color = [0 0 0];

% Offset error bars slightly from bar centers for better visibility
er(1).XData = er(1).XData - 0.15; % b(1).BarWidth * 0.5;
er(2).XData = er(2).XData + 0.15; %;b(1).BarWidth * 0.5;

% Label x-axis with ROI names at 45 degree angle
xticks(1:length(roiname))
xticklabels(roiname)
xtickangle(45)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate average values and profiles for each ROI and contrast
% This section processes and organizes brain imaging data by:
% - Generating group-level and age-grouped profiles for different contrasts
% - Calculating average projections and skewness measures
% - Processing global brain profiles
% - Preparing matrices for structural measurements (thickness, surface, volume)
% All results are saved in CSV files and organized data tables for further analysis


%% Calculate average values and profiles for each ROI and contrast
% Create output directory for group results
system(['mkdir ', outdir_group]);
depth_sampled = 0.05:0.05:0.95;
close all

% Generate profile files and plots
% Make group profiles for selected contrasts
function_to_make_group_profiles(freesurfer_folder, contrasts([1:3,5:7]), depth_sampled, atlas, nodenr, outdir_ind, outdir_group)

% Generate age-grouped profile plots
function_to_make_age_grouped_profiles(freesurfer_folder, contrasts([1:3,5:7]), depth_sampled, atlas, nodenr, outdir_ind, outdir_group)

% Calculate average projections and load values
[average_projection average_skewness all_table all_table_iqr all_table_skewness all_table_skewness_paquola] = function_make_average_projection(contrasts([1:3,5:7]), atlas, nodenr, outdir_group);

% Merge skewness measures into main table
fn = fieldnames(all_table_skewness);
for v = 1:length(fn)
   % Add regular skewness measure
   newname = ['skewness_',fn{v}];
   all_table.(newname) = all_table_skewness.(fn{v});
   
   % Add Paquola skewness measure
   newname = ['skewness_paquola_',fn{v}];
   all_table.(newname) = all_table_skewness_paquola.(fn{v});
end

% Commented out test code for skewness correlation plots
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Analyze and visualize relationships between whole brain measures
% Creates a figure with 3 scatter plots showing relationships between:
% 1. Brain weight vs intracranial volume (ICV)
% 2. Brain weight vs cortical volume 
% 3. White matter volume vs cortical volume
% Each plot includes a least squares regression line and data point labels

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Model and plot whole brain measures across age
% Creates a figure with subplots for each whole brain measure (from wb struct)
% Each subplot shows:
% - Scatter plot of the measure vs age
% - Fitted Hallgren growth curve (exponential model)
% Saves figure as 'Lineplot_wholebrain_measures.png'

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ROI-based Analysis
% This section performs region of interest (ROI) based analysis on brain imaging data
% Key operations:
% - Processes data for baby scans separately and stores median values
% - Handles multiple contrasts (R1, MTsat, etc.) for each ROI
% - Calculates median values and interquartile ranges
% - Saves results in tables for further analysis

% Define ROIs to analyze (first 38 regions)
mask_idx = 1:38;

% Process baby scans if not already done
babyfile = [resultdir,'/BabyROIvals_median.mat'];
if exist(babyfile) == 0
    % Initialize counter for baby scans
    count2 = 1;
    
    % Process specific baby scans (indices 2 and 9)
    for scan2 = [2 9]
        % Get brain ID and format for filenames
        brain_slash_id = brain_slash_ids{scan2};
        brain_undersc_id = strrep(brain_slash_id,'/','_'); 

        % Get age information
        id = brain_slash_ids{scan2};
        idx = str2num(id(1:3));
        ages2(count2) = ages_from_database(idx);

        % Define paths for ROI and MPM data
        ROI_folder = ['preprocessed/',brain_slash_id,'_V2/manual_ROIs'];
        MPM_folder = ['preprocessed/',brain_slash_id,'_V2/MPMs_to_use'];
        
        % Process each contrast
        for c = 1:4
            contrast = contrasts{c};
            % Get MPM file path
            file = [MPM_folder,'/',contrast,'_0p3_run01_brain_masked_reoriented.nii']
            
            % Calculate median and IQR statistics using FSL
            [error val] = system(['fslstats ', file, ' -P 50']);  % Median
            [error valu] = system(['fslstats ', file, ' -P 75']); % 75th percentile
            [error vall] = system(['fslstats ', file, ' -P 25']); % 25th percentile
            
            % Store values in tables
            all_baby_table.(contrast)(:,count2) = repmat(str2double(val),nodenr*2,1);
            all_baby_iqr_table.(contrast)(:,count2) = repmat(str2double(valu)-str2double(vall),nodenr*2,1);
        end
        count2 = count2 + 1;
    end
    % Save baby scan results
    save(babyfile,'all_baby_table','all_baby_iqr_table','ages2');
else
    % Load existing baby scan results
    load(babyfile,'all_baby_table','all_baby_iqr_table','ages2');
end

% Get list of all variables for processing
all_variables = fieldnames(all_table);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Process and organize brain imaging data by contrast and region
% - Extracts values for left and right hemispheres
% - Incorporates baby scan data if available
% - Calculates medians and interquartile ranges (IQR)
% - Sorts data by age
% - Computes means across hemispheres

cs = 1:length(all_variables);
value_name = 'test';
clear vals vals_left vals_right vals_mean vals_iqr_left vals_iqr_right vals_iqr_mean

for c = cs
    contrast = all_variables{c};
       for region = 1:length(mask_idx)
            % Get original values for this contrast
            origvals = all_table.(contrast); 
            
            % Process median values
            try
                % Add baby scan data if available
                vals_baby = all_baby_table.(contrast)';
                origvals = [origvals; vals_baby];
                [ages_sorted, sort_idx] = sort([ages, ages2]);
            catch
                % Use only main dataset if no baby data
                [ages_sorted, sort_idx] = sort([ages]); 
            end
            
            % Store left and right hemisphere values
            vals.left.(contrast)(:,region) = origvals(sort_idx,mask_idx(region));
            vals.right.(contrast)(:,region) = origvals(sort_idx,mask_idx(region)+38);
            
            % Process IQR values
            try
                % Add baby scan IQR data if available
                vals_iqr = all_table_iqr.(contrast); 
                vals_baby_iqr = all_baby_iqr_table.(contrast)';
                vals_iqr = [vals_iqr; vals_baby_iqr];
                [ages_sorted, sort_idx] = sort([ages, ages2]);
            catch
                % Use zeros for IQR if no data available
                vals_iqr = zeros(size(origvals));
                [ages_sorted, sort_idx] = sort([ages]);
            end
            
            % Store left and right hemisphere IQR values
            vals.iqr_left.(contrast)(:,region) = vals_iqr(sort_idx,mask_idx(region));
            vals.iqr_right.(contrast)(:,region) = vals_iqr(sort_idx,mask_idx(region)+38);
            vals.age.(contrast) = ages_sorted';
       end
       
       % Calculate means across hemispheres
       vals.mean.(contrast) = nanmean(cat(3 ,vals.left.(contrast), vals.right.(contrast)), 3); 
       vals.iqr_mean.(contrast) = nanmean(cat(3 ,vals.iqr_left.(contrast), vals.iqr_right.(contrast)), 3); 
end        

%% Fit models to the data
% This section fits two types of models to the data:
% 1. A Hallgren-style exponential model for most contrasts
% 2. A linear model for surface, volume and skewness contrasts
% Models are fit separately for left, right and mean values across hemispheres
% Parameters and fit statistics are stored in model and model2 structs

% Initialize model structures
clear model model2
hemispheres = {'left','right','mean'}

% Loop through brain regions, hemispheres and contrasts
for region = 1:38 %length(mask_idx)
    for h = 1:length(hemispheres)
        hemi = hemispheres{h};
         for c = cs 
            contrast = all_variables{c};
            values_to_use = vals.(hemi).(contrast)(:,region);
            
            % Set up variables for model fitting
            iv = vals.age.(contrast); % Independent variable (age)
            dv = values_to_use; % Dependent variable (contrast values)
            fake_iv = 0:0.01:max(iv); % Generate dense x-axis for predictions
            
            % Define Hallgren exponential model
            newfun = @(cpars,iv)(cpars(1)*(1-exp(cpars(2)*iv*-1))+cpars(3));
            pars0 = [5, 0.1, 1]; % Initial parameters: amplitude, rate, constant; estimated from paper
            
            % Use modified exponential decay for thickness
            if strcmp(contrast,'thickness') 
                newfun = @(cpars,iv)(cpars(1)*(exp(cpars(2)*iv*-1))+cpars(3));
            end
            warning('') % Clear previous warnings
            %% fit linear model
            if strcmp(contrast,'surface') || strcmp(contrast,'volume') || strcmp(contrast(1:2),'sk')   %%% model with a linear model 

                pars0 = [0 0];
                newfun = @(cpars,iv)(cpars(2) + cpars(1)*iv); % Linear function
                
                % Remove NaN values
                nonnanidx = find(~isnan(dv));
                
                % Set parameter bounds based on contrast type
                if strcmp(contrast,'surface') || strcmp(contrast,'volume')
                   lb = [-9999, 0]; ub = [9999 9999]; % Constrain constant to be positive
                elseif strcmp(contrast(1:2),'sk')
                   lb = [-9999, -9999]; ub = [9999 9999]; % No constraints for skewness
                end 
                
                % Fit model using least squares
                [BETA,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(newfun,pars0,iv(nonnanidx),dv(nonnanidx),lb,ub);
                
                % Calculate confidence intervals and predictions
                ci = nlparci(BETA, RESIDUAL,'Jacobian',JACOBIAN);
                prediction = newfun(BETA, iv);
                fakeprediction = newfun(BETA, fake_iv);
                
                % Add placeholder amplitude parameter
                BETA = [0 BETA];
                ci = [0 0; ci];
                
            else
                % Fit exponential model
                try
                    nonnanidx = find(~isnan(dv));
                    
                    % Set parameter bounds based on contrast type
                    if strcmp(contrast(1:2),'th')
                        lb = [-9999, 0, 0]; ub = [5, 9999, 5]; % Constrain thickness < 5mm
                    elseif strcmp(contrast(1:2),'sk')
                        lb = [-10, -9999, -10]; ub = [10, 9999, 10]; % Bounds for skewness
                    else
                        lb = [0, 0, 0]; ub = [9999, 9999, 9999]; % Non-negative bounds for MPMs
                    end
                    
                    % Fit model using least squares
                    [BETA,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = lsqcurvefit(newfun,pars0,iv(nonnanidx),dv(nonnanidx),lb,ub);
                    
                    % Calculate confidence intervals and predictions
                    ci = nlparci(BETA, RESIDUAL,'Jacobian',JACOBIAN);
                    prediction = newfun(BETA, iv);
                    fakeprediction = newfun(BETA, fake_iv);
                catch
                    % Handle fitting errors
                    BETA = NaN(1,length(pars0));
                    ci = NaN(length(pars0),2);
                    prediction = NaN(size(iv));
                    fakeprediction = NaN(size(fake_iv));
                end
            end
           
            % Store confidence intervals
            model.(contrast).(hemi).rate_CI_low(region) = ci(2,1);
            model.(contrast).(hemi) = model.(contrast).(hemi);
            model.(contrast).(hemi).rate_CI_high(region) = ci(2,2);    
            model.(contrast).(hemi).amplitude_CI_low(region) = ci(1,1);
            model.(contrast).(hemi).amplitude_CI_high(region) = ci(1,2);       
            model.(contrast).(hemi).constant_CI_low(region) = ci(3,1);
            model.(contrast).(hemi).constant_CI_high(region) = ci(3,2); 

            % Calculate fit statistics
            [r p] = corrcoef(dv,prediction,'rows','complete');
            model.(contrast).(hemi).coeffr2(region) = r(1,2).^2;
            model.(contrast).(hemi).p(region) = p(1,2);
            model.(contrast).(hemi).n(region) = sum(~isnan(dv));
            
            % Store model parameters
            model.(contrast).(hemi).rate(region) = BETA(2);
            model.(contrast).(hemi).amplitude(region) = BETA(1);
            model.(contrast).(hemi).constant(region) = BETA(3);

            % Check for warnings during fitting
            [warnMsg, warnId] = lastwarn;
            if ~isempty(warnMsg)
                1 + 1 
                model.(contrast).(hemi).problem(region) = 1;
            else
                model.(contrast).(hemi).problem(region) = 0;
            end
            
            % Store raw data and predictions
            model.(contrast).(hemi).iv(:,region) = iv;
            model.(contrast).(hemi).dv(:,region) = dv;
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


%% Plot trajectories for selected brain regions
% This section creates line plots comparing model fits across 6 key brain regions:
% Visual, Auditory, Motor, Somatosensory, Broca homologue, and Prefrontal areas.
% For each contrast, it generates two figures:
% 1) Hallgren model fits with exponential/linear trajectories
% 2) Miller model fits with quadratic trajectories
% The plots include raw data points and fitted curves with legends showing model parameters.

%% Define regions and colors
% Define names and indices for 6 key brain regions to analyze
mask_names2 = {'Visual','Auditory','Motor','Somatosensory','Broca homologue','Prefrontal'};
mask_idx2 = [24 34 1 26 5 8]; 
% Define color codes for plotting each region
six_color_codes = {[140 20 60], [26 58 148], [60 140 180], [20 180 140], [220 20 100], [180 22 140]};


%% Generate line plots
hemi = 'mean'
% Loop through selected contrast indices
for c = [1:3, 7, 9, 11, 19] %1:3 %cs %s = 1:length(all_variables)
    % Create two figures - one for Hallgren model, one for Miller model
    f1 = figure();
    set(f1,'Position',[100 100 1900 700]);
    f2 = figure(); 
    set(f2,'Position',[100 100 1500 700]);
    
    % Loop through each brain region
    for regiontoplot = 1:length(mask_idx2)
        contrast = all_variables{c};
        region_idx = mask_idx2(regiontoplot);
        
        % Get data for current region from Hallgren model
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
        
        %% Plot Hallgren model results
        figure(f1)
        hold on
        % Plot raw data points
        scatter(iv, dv, 100, six_color_codes{regiontoplot}./255,'filled');
        
        % Create legend entry based on model type and significance
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
        
        % Plot fitted curve
        aah = plot(fake_iv,fake_pred,'Color',six_color_codes{regiontoplot}./255,'LineWidth',3.5);
        aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold off
        
        % Set axis labels and formatting
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
        
        % Set y-axis limits if all values are positive
        all_min(regiontoplot) = min([dv]);
        all_max(regiontoplot) = max([dv]);
        if min(all_min) > 0 %%% only if no negative values
            try
                m = max(all_max);
                ylim([0, max(m) + 0.2*max(m)]);           
            end
        end
        
        %% Plot Miller model results
        % Get data for current region from Miller model
        thismodel2 = model2.(contrast).(hemi);
        iv = thismodel2.iv(:,region_idx);
        dv = thismodel2.dv(:,region_idx);
        pred = thismodel2.pred(:,region_idx);
        fakepred = thismodel2.fakepred(:,region_idx);
        fake_iv = thismodel2.fake_iv(:,region_idx);
        r2val = thismodel2.r2(:,region_idx);
        
        figure(f2)
            hold on;
            % Plot raw data and fitted curve
            scatter(iv,dv, 100, six_color_codes{regiontoplot}./255,'filled');
            h = plot(fake_iv,fakepred,'Color',six_color_codes{regiontoplot}./255,'LineWidth',3.5);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            hold off
            
            % Set axis labels and formatting
            xlabel('years');
            ylabel(contrast, 'Interpreter', 'none');
            set(gca, 'FontSize',25);
            set(gca, 'Color', 'white');
            legend_entry2{regiontoplot} = [mask_names2{regiontoplot},': ', num2str(round(r2val,2))]; 
        
    end
 
    % Save figures with legends
    figure(f1)
        legend(legend_entry,'Location','NorthWestOutside','FontSize',15);
        saveas(gcf,[resultdir,'/figures/Lineplot_',contrast,'_',hemi,'.png']);
    figure(f2)
        legend(legend_entry2,'Location','NorthWestOutside','FontSize',15);
        saveas(gcf,[resultdir,'/figures/Lineplot_QuadraticFit_',contrast,'_',hemi,'.png']);  

end




%% Create individual trajectory plots for each region and contrast
% This section creates subplots showing age-related trajectories for each brain region
% separately for R1, R2* and MTsat contrasts. Each subplot shows raw data points
% and fitted exponential curves.

hemi = 'mean'
for c = 1:3 
    contrast = all_variables{c};
    % Create new figure with specified dimensions
    f1 = figure();
    set(f1,'Position',[100 100 1900 700]);
    subplcount = 1;
    
    % Loop through all 38 regions
    for region_idx = 1:38
        % Get model data for current region and contrast
        thismodel = model.(contrast).(hemi);
        iv = thismodel.iv(:,region_idx);          % Independent variable (age)
        dv = thismodel.dv(:,region_idx);          % Dependent variable (contrast value) 
        fake_iv = thismodel.fake_iv(:,region_idx); % Interpolated x values for smooth curve
        fake_pred = thismodel.fake_pred(:,region_idx); % Predicted values for smooth curve
        pred = thismodel.pred(:,region_idx);
        rate = thismodel.rate(:,region_idx);
        constant = thismodel.constant(:,region_idx);
        amplitude = thismodel.amplitude(:,region_idx);
        pval = thismodel.p(:,region_idx);
        r2 = thismodel.coeffr2(:,region_idx);
        
        % Create subplot for current region
        figure(f1)
        subplot(6,7,subplcount); subplcount = subplcount + 1;
        hold on
        scatter(iv, dv, 30, six_color_codes{regiontoplot}./255,'filled');
        aah = plot(fake_iv,fake_pred,'Color',six_color_codes{regiontoplot}./255,'LineWidth',1.5);
        aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold off
        xlabel('years');
        if strcmp(contrast,'R1')
            conlab = 'R1 (s^{-1})';
        elseif strcmp(contrast,'R2s')
            conlab = 'R2* (s^{-1})';
        elseif strcmp(contrast,'MTsat')
            conlab = 'MTsat (p.u.)';
        else
            conlab = contrast;
        end
        ylabel(conlab);
        
        % Set plot formatting
        set(gca, 'FontSize',10);
        set(gca, 'Color', 'white');
        
        % Adjust y-axis limits for positive values
        all_min(regiontoplot) = min([dv]);
        all_max(regiontoplot) = max([dv]);
        if min(all_min) > 0 % Only adjust if no negative values
            try
                m = max(all_max);
                ylim([0, max(m) + 0.2*max(m)]);           
            end
        end        
    end
end

%% Create output tables with model results
% This section creates tables containing model fit parameters and statistics for either:
% 1) 6 key ROIs (Visual, Auditory, Motor, Somatosensory, Broca, Prefrontal)
% 2) All 38 cortical regions
% Tables include sample size, R^2, p-values, model parameters and confidence intervals

% Loop through table types (6 ROIs or all regions)
table_types = {'6rois','all'};
for tabtodo = 1:2
    % Set up region indices and names based on table type
    if strcmp(table_types{tabtodo}, 'all')
        % For all regions table, use full 38 region set
        mask_idx2 = 1:38;
        mask_names2 = num2cell([1:38]);
        % Get region names from annotation file
        annotfile = ['results/Cortical_analysis/group_freesurfer/fsaverage/label/lh.BB38chimp.annot'];
        [vertices, labeling, colortable] = read_annotation(annotfile);
        mask_names2 = colortable.struct_names(2:end)';
    else
        % For 6 ROIs table, use key regions
        mask_names2 = {'Visual','Auditory','Motor','Somatosensory','Broca homologue','Prefrontal'};
        mask_idx2 = [24 34 1 26 5 8]; 
    end
    
    % Set up table structure
    firstcol = mask_names2';
    mat1 = [];
    
    % Loop through contrasts and compile results
    for c = [1:3, 7, 9, 11, 19] %cs
        contrast = all_variables{c};
        shortcut = model.(contrast).(hemi);
        miller_r2 = model2.(contrast).(hemi).r2;
        
        % Calculate percent of final value present at birth
        perc_at_birth = 100 * shortcut.constant(mask_idx2) ./ (shortcut.constant(mask_idx2) + shortcut.amplitude(mask_idx2));
        
        % Compile all metrics into matrix
        mat1 = [shortcut.n(mask_idx2)', ...                    % Sample size
               round(shortcut.coeffr2(mask_idx2),2)', ...      % R-squared
               round(shortcut.p(mask_idx2),3)', ...            % p-value
               round(shortcut.constant(mask_idx2),2)', ...     % Constant term
               round(shortcut.amplitude(mask_idx2),2)', ...    % Amplitude
               round(shortcut.constant(mask_idx2) + shortcut.amplitude(mask_idx2),2)', ... % Plateau
               round(perc_at_birth)', ...                      % Percent at birth
               round(1./shortcut.rate(mask_idx2),2)', ...      % Time constant
               round(shortcut.rate(mask_idx2),4)', ...         % Rate
               round(shortcut.rate_CI_low(mask_idx2),2)', ...  % Lower CI
               round(shortcut.rate_CI_high(mask_idx2),2)', ... % Upper CI
               round(miller_r2(mask_idx2)',2)];                % Miller model R2
               
        % Create and write table
        table_to_write = cell2table([firstcol,num2cell(mat1)],'VariableNames',...
            {'Region','n','R^2', 'p', 'constant', 'amplitude', 'plateau', 'perc', ...
             'time const', 'rate', 'lower rate CI', 'higher rate CI','MillerR2'});
        writetable(table_to_write,[resultdir,'/model_results_',contrast,'_',hemi,'_',table_types{tabtodo}]);
    end  
end



%%% average time constant values across cortex for mtsag
median(1./model.MTsat.mean.rate)
3.1028 * 3 %%% highest realistic value x 3


%%% Analyze and plot relationships between rate/time constant and amplitude
% This section examines correlations between the rate of development (and time constant)
% and the final amplitude across brain regions for R2* and MTsat measures.
% For each measure, it creates scatter plots showing:
% 1) Rate vs amplitude 
% 2) Time constant vs amplitude
% Only includes regions with significant model fits (p<0.05) and removes outliers.

varstotake = {'R2s','MTsat'};
for v = 1:2
    f = figure()
    set(f,'Position',[100 100 300 300])
    var = varstotake{v};
    
    % Get model parameters for left and right hemispheres
    pvals = [model.(var).left.p, model.(var).right.p];
    a1 = [1./model.(var).left.rate, 1./model.(var).right.rate]; % Time constants
    a2 = [model.(var).left.rate, model.(var).right.rate]; % Rates
    b = [model.(var).left.constant + model.(var).left.amplitude, model.(var).right.constant + model.(var).right.amplitude]; % Final amplitudes
    
    % Plot Rate vs Amplitude
    idx = find(pvals < .05 & ~isoutlier(a1) & ~isoutlier(b)); % Select significant regions, remove outliers
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
    
    % Plot Time Constant vs Amplitude
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


%% Assess parameter reliability by comparing left vs right hemispheres
% This section examines the reliability of model parameters by correlating rates
% between left and right hemispheres for each measure (MTsat, R2s, surface, volume, 
% thickness, skewness_R2s). Creates a 3x2 subplot figure showing scatter plots
% with fitted lines for each measure. Strong correlations indicate reliable parameters.

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

%% Analyze correlations between microstructural parameters across regions
% Create data structures containing key model parameters for each measure:
% - Time constants (1/rate)
% - Rates 
% - Final amplitudes (constant + amplitude)
% Also includes skewness measures and surface expansion metrics
variablestruc = {[1./model.R2s.left.rate, 1./model.R2s.right.rate], [model.R2s.left.rate, model.R2s.right.rate],[model.R2s.left.constant + model.R2s.left.amplitude, model.R2s.right.constant + model.R2s.right.amplitude] ...
    [1./model.MTsat.left.rate, 1./model.MTsat.right.rate],[model.MTsat.left.rate, model.MTsat.right.rate],  [model.MTsat.left.constant + model.MTsat.left.amplitude, model.MTsat.right.constant + model.MTsat.right.amplitude] ...
    median(all_table_skewness.R2s), median(all_table_skewness.MTsat),nanmedian(all_table_skewness.R1),  ontogenetic_surface_expansion_perc}; %, phylogenetic_surface_expansion};

% Define variable names corresponding to each parameter
variable_names = {'R2s_time_constant','R2s_rate','R2s_amplitude','MTsat_time_constant','MTsat_rate','MTsat_amplitude','R2s_skewness','MTsat_skewness','R1_skewness','ontogenetic_surface_expansion'}; %,'phylogenetic_surface_expansion' };

%% FreeSurfer annotation commands for reference
% Commands to run in terminal to generate template annotations:
% ['mris_ca_label -t ',atlas_folder,'/BB38chimp.annot.ctab fsaverage lh ',freesurfer_folder,'/surf/lh.sphere.reg ', atlas_folder,'/lh.BB38chimp.gcs ',freesurfer_folder,'/label/lh.BB38chimp.annot']
% ['mris_ca_label -t ',atlas_folder,'/BB38chimp.annot.ctab fsaverage rh ',freesurfer_folder,'/surf/rh.sphere.reg ', atlas_folder,'/rh.BB38chimp.gcs ',freesurfer_folder,'/label/rh.BB38chimp.annot']

%% Generate surface projections for each variable
% Loop through variables and create surface plots showing parameter distributions
for var = 1:length(variable_names)
    
    variable_name = variable_names{var};
    variable_vec = variablestruc{var};
    figure_file_to_save_prefix = [resultdir,'/roi_wise_projection_figures/fsaverage_',hemisphere,'_',variable_name,'_brewermap'];
    projection_file_to_save = [resultdir,'/roi_wise_projection_figures/fsaverage_',hemisphere,'_',variable_name,'.mgh'];
    
    function_make_roiwise_surface_plot(variable_vec, variable_name, [1 99], projection_file_to_save, figure_file_to_save_prefix)

end

% Clean up and save workspace
close all
save([resultdir,'/microstructure_pipeline_workspace'])
%adult_surf_surface_matrix(:,ages>17)

