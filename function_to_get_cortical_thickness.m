function [thickness, surface, volume, normalised_volume, roiname, global_vals] = function_to_get_cortical_thickness(freesurfer_folder)

    thickness = []; surface = []; volume = []; ICV = []; normalized_cortex_gm_volume = [];
    hemispheres = {'lh','rh'};
    for hem = 1:length(hemispheres)
            clear C
            BB38_stats_file = [freesurfer_folder,'/stats/',hemispheres{hem},'.BB38chimp.stats'];
            if exist(BB38_stats_file) == 2
                fid = fopen(BB38_stats_file);
                C = textscan(fid,'%s %f %f %f %f %f %f %f %f %f', 'headerlines', 60);
                thickness = [thickness; C{5}]; %%% in mm
                surface = [surface; C{3}]; %%% in mm2
                volume = [volume; C{4}]; %%% in mm3
                roiname = C{1};
            end
    end
    
    %statsfile = [freesurfer_folder,'/stats/',hemispheres{hem},'.aparc.stats'];

    statsfile = [freesurfer_folder,'/stats/brainvol.stats'];

if exist(statsfile) == 2

    %fid = fopen(statsfile);
    %clear C
    %C = textscan(fid, '%f %f %f %f %s %f %f %f %f %f', 'headerlines', 79);
    %fclose(fid);
    %%% these need to be corrected for intracranial.. while
    %%% cortical thickness does not need to be: https://surfer.nmr.mgh.harvard.edu/fswiki/eTIV
    %[a b] = system(['cat ',statsfile,' | grep "Measure EstimatedTotalIntraCranialVol, eTIV, Estimated Total Intracranial Volume"']);
    [a b] = system(['cat ',statsfile,' | grep "Measure BrainSeg, BrainSegVol, Brain Segmentation Volume"']);
    ss1 =  strsplit(b,'ume,');
    ss2 = strsplit(ss1{2}, ', mm');
    global_vals.ICV = str2double(ss2{1}); %[ICV; str2num(b(85:length(b)-7))];
    
    [a b] = system(['cat ',statsfile,' | grep "Measure Cortex, CortexVol, Total cortical gray matter volume"']);
    ss1 =  strsplit(b,'ume,');
    ss2 = strsplit(ss1{2}, ', mm');
    global_vals.cortical_vol = str2double(ss2{1}); 
    
    normalised_volume = 100 * volume / global_vals.cortical_vol; %%% in percent from cortical volume
    
    [a b] = system(['cat ',statsfile,' | grep "Measure CerebralWhiteMatter, CerebralWhiteMatterVol, Total cerebral white matter volume"']);
    ss1 =  strsplit(b,'ume,');
    ss2 = strsplit(ss1{2}, ', mm');
    global_vals.wm_vol = str2double(ss2{1});  
    
    %[a b] = system(['cat ',statsfile,' | grep "Measure Cortex, CortexVol, Total cortical gray matter volume"']);
    %normalized_cortex_gm_volume = [normalized_cortex_gm_volume; 100 * (str2num(b(64:length(b)-7)) / ICV)];
   % normalized_cortex_gm_volume = NaN;
end