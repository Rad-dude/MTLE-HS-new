% RKUL_MTLE-HS-preproc_new v 0.1 - dd 3/11/2018 @ radwan - sunaert We will
% be doing both the T1 and DTI preproc here

%%%% We need to split our atlases into 4D and use em like that, it's
%%%% cleaner, much cleaner...

%% Part 1 Configure my stuff...
clear all;
clc;
dir_main   = '/Volumes/LaCie/MTLE_HS';
ncpu = '7';

% define your dirs & make them
dir_source_DTI = [dir_main filesep 'NII' filesep 'DTI' filesep 'orig'];
dir_source_T1 = [dir_main filesep 'NII' filesep 'T1s' filesep 'orig'];
dir_t1_int_out = [dir_main filesep 'NII' filesep 'T1s' filesep 'intermediate'];
dir_t1_out = [dir_main filesep 'NII' filesep 'T1s' filesep 'T1_output'];
% dir_t1_temp_wskull = [dir_t1_int_out filesep 'T1_template_wskull'];
% dir_t1_temp_brain = [dir_t1_int_out filesep 'T1_template_brain'];
dir_DTI_int_out = [dir_main filesep 'NII' filesep 'DTI' filesep 'intermediate'];
dir_DTI_out = [dir_main filesep 'NII' filesep 'DTI' filesep 'DTI_output'];
dir_log    = [dir_main filesep 'LOG'];
dir_templates = [dir_main filesep 'templates'];
mkdir(dir_log);
mkdir(dir_t1_int_out);
mkdir(dir_t1_out);
mkdir(dir_DTI_int_out);
mkdir(dir_DTI_out);

% list the T1s found and the DTIs found
diary([dir_log filesep 'command_log_dwi_preproc_2.txt']);
tic
s_DTI = (dir([dir_source_DTI filesep '*.nii.gz']));
s_T1s = (dir([dir_source_T1 filesep '*_orig.nii.gz']));
DTI_files = (extractfield(s_DTI,'name'))';
T1s_files = (extractfield(s_T1s,'name'))';

unix(['touch ' dir_log filesep 'sanity_check_fail.txt']);
unix(['touch ' dir_log filesep 'sanity_check_pass.txt']);

% quick sanity check to see if all T1s are paired with a DTI
for i = 1:(size(T1s_files,1))
    crs_check = strncmpi(T1s_files(i), DTI_files(i), 6);
    if crs_check == 0
        unix(['echo "*** It seems that ' char(T1s_files(i)) ' has no DWI data" >> ' dir_log filesep 'sanity_check_fail.txt' ]);
    else
        unix(['echo "*** It seems that ' char(T1s_files(i)) ' has DWI data" >> ' dir_log filesep 'sanity_check_pass.txt' ]);
    end
end

% Defining template images and other variables
T1_template = [dir_templates filesep 'PD25-T1MPRAGE-template-1mm.nii.gz'];
T1_template_brain = [dir_templates filesep 'PD25-T1MPRAGE-template-1mm_brain.nii.gz'];
T1_temp_brain_mask = [dir_templates filesep 'PD25-atlas-mask-1mm.nii.gz'];
T2_template = [dir_templates filesep 'mni_icbm152_t2_tal_nlin_asym_09a.nii'];
T2_temp_brain_mask = [dir_templates filesep 'mni_icbm152_t2_tal_nlin_asym_09a_mask.nii'];
T2_temp_brain = [dir_templates filesep 'mni_icbm152_t2_tal_nlin_asym_09a_brain.nii.gz'];
aal_template_brain = [dir_templates filesep 'MNI152_T1_1mm_brain.nii.gz'];
aal_labels = [dir_templates filesep 'aal_atlas.nii'];
MNI_lobes_rl = [dir_templates filesep 'MNI_rl1t.nii.gz'];
% The manually defined exclude ROIs exclude_labels = [dir_templates filesep
% 'exclude_ROIs.nii.gz']; txt_aal_labels = [dir_templates filesep
% 'aal_atlas.txt'];
HOCA = '$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr25-1mm.nii.gz';
HOSCA = '$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-1mm.nii.gz';
JHU_WM = '$FSLDIR/data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr0-1mm.nii.gz';
s_dg_DTI = (dir([dir_DTI_int_out filesep '*' filesep '*_dn_dg.mif']));
s_pp_DTI = (dir([dir_DTI_int_out filesep '*' filesep '*_pp.mif']));


%% Part 2 prep the T1s and DTIs up to mrdegibbs

for i = 1:(size(T1s_files,1))
    % defining all variables for the first loop - this goes up to degibbs
    T1_subj = char(T1s_files(i));
    DTI_subj = char(DTI_files(i));
    subj_basename = (T1_subj(1:end-18));
    T1_subj_reori = ([subj_basename 'T1_reori.nii.gz']);
    T1_subj_1mm = ([subj_basename 'T1_1mm.nii.gz']);
    T1_subj_bc = ([subj_basename 'T1_bc.nii.gz']);
    T1_subj_brain = ([subj_basename 'T1_BrainExtractionBrain.nii.gz']);
    subj_t1_int_out = [dir_t1_int_out filesep subj_basename 'T1_int_out'];
    subj_t1_out  = [dir_t1_out filesep subj_basename 'T1_output'];
    subj_DTI_int_out = [dir_DTI_int_out filesep subj_basename 'DTI_int_out'];
    %     subj_DTI_out  = [dir_DTI_out filesep subj_basename 'DTI_output'];
    DTI_sub_base = char(DTI_files(i));
    DTI_sub_base = (DTI_sub_base(1:end-7));
    DTI_bval = [DTI_sub_base '.bval']; % diff weighting input
    DTI_bvec = [DTI_sub_base '.bvec'];
    DTI_bmatrix = [DTI_sub_base '_bmatrix.txt']; % diff weighting int.
    DTI_sub_nii = [DTI_sub_base '.nii.gz'];
    DTI_sub_mif = [DTI_sub_base '.mif'];
    DTI_1 = [DTI_sub_base '_dn.mif'];
    DTI_2 = [DTI_sub_base '_dn_dg.mif'];
    
    subj_dti_int_out = [dir_DTI_int_out filesep subj_basename 'DTI_int_out'];
    subj_dti_out = [dir_DTI_out filesep subj_basename 'DTI_output'];
    mkdir(subj_t1_int_out);
    mkdir(subj_t1_out);
    mkdir(subj_dti_int_out);
    mkdir(subj_dti_out);
    
    % First we reorient the image to standard orientation with FSL
    unix(['source ~/.bash_profile ; fslreorient2std ' dir_source_T1 filesep T1_subj ' ' subj_t1_int_out filesep T1_subj_reori]);
    
    % Second we use mrresize to put the T1s to 1 mm isotropic dimensions
    unix(['source ~/.bash_profile ; mrresize -force -nthreads 7 -voxel 1 ' subj_t1_int_out filesep T1_subj_reori ' ' subj_t1_int_out filesep T1_subj_1mm]);
    
    % ANTS N4 bias field correction
    unix(['source ~/.bash_profile ; N4BiasFieldCorrection -d 3 -i ' subj_t1_int_out filesep T1_subj_1mm ' -o ' subj_t1_int_out filesep T1_subj_bc]);
    
    % ANTS BET % there was a typo here resulting in names like this
    % PT_001__T1_, this is now solved.
    unix(['source ~/.bash_profile ; antsBrainExtraction.sh -d 3 -a ' subj_t1_int_out filesep T1_subj_bc ' -e ' T1_template ' -m ' T1_temp_brain_mask ' -o ' subj_t1_int_out filesep subj_basename 'T1_ -s nii.gz -u 1']);
    
    % Convert DTI nii.gz, .bvals and .bvecs to .mif and .txt
    unix(['source ~/.bash_profile ; mrconvert -nthreads 7 -force -fslgrad ' dir_source_DTI filesep DTI_bvec ' ' dir_source_DTI filesep DTI_bval ...
        ' ' dir_source_DTI filesep DTI_sub_nii ' ' subj_DTI_int_out filesep DTI_sub_mif ' -export_grad_mrtrix ' subj_DTI_int_out filesep DTI_bmatrix]);
    
    % Denoise and degibbs
    unix(['source ~/.bash_profile ; dwidenoise -force -nthreads 7 ' subj_DTI_int_out filesep DTI_sub_mif ' ' subj_DTI_int_out filesep DTI_1]); % step 1 denoise
    
    unix(['source ~/.bash_profile ; mrdegibbs -force -nthreads 7 ' subj_DTI_int_out filesep DTI_1 ' ' subj_DTI_int_out filesep DTI_2]);
end

%% Part 3 DTI preprocessing with parfor and Eddy_cuda, then bias correct in a regular for loop



poolobj = gcp('nocreate');
delete(poolobj)
parpool(4);
parfor i = 1:(size(s_dg_DTI,1))
    DTI_3 = [s_dg_DTI(i).name(1:end-9) 'pp.mif'];
    
    unix(['source ~/.bash_profile ; time dwipreproc -force -nthreads 2 -rpe_none -pe_dir AP -eddy_options " --repol --residuals" ' ...
        s_dg_DTI(i).folder filesep s_dg_DTI(i).name ' ' s_dg_DTI(i).folder filesep DTI_3]); % added time ahead of dwipreproc for timing this step
    
end

for i = 1:(size(s_pp_DTI,1))
    DTI_4 = [s_pp_DTI(i).name(1:end-6) 'pp_bc.mif'];
    
    unix(['source ~/.bash_profile ; dwibiascorrect -force -ants -nthreads 7 ' s_pp_DTI(i).folder filesep s_pp_DTI(i).name ' ' s_pp_DTI(i).folder filesep DTI_4]);
    
end

%% Part 4 DWI tensor, maps, response function, SH and peaks calculation, and VOIs prep

s_bc_DTI = (dir([dir_DTI_int_out filesep '*' filesep '*_pp_bc.mif']));
DTI_outs  = (dir([dir_DTI_out filesep 'PT_*_DTI_output']));

% Define VOIs to use from whatever atlas
VOIs = struct('name', {'056 - PCC', '057 - PC', '058 - R_Cuneal', '059 - L_Cuneal', '062 - R_aPaHC', ...
    '063 - L_aPaHC', '064 - R_pPaHC', '065 - L_pPaHC', ...
    '100 - R_Hipp', '101 - L_Hipp', '102 - R_Amyg', '103 - L_Amyg', '013 - Rt_PreCG', '014 - Lt_PreCG',  ... 
    '001 - exclude_Lt_WM', '008 - BStem', '012 - exclude_Rt_WM', '002 - exclude_Lt_GM', '013 - exclude_Rt_GM', ...
    '008 - rt_temporal', '080 - lt_temporal', '037 - Rt_VCing', '038 - Lt_VCing'}, 'lt', {'55.5', '56.5', '57.5', '58.5', ...
    '61.5', '62.5', '63.5', '64.5', '99.5', '100.5', '101.5', '102.5', '12.5', '13.5', '0.5', '7.5', '11.5', '1.5', ...
    '12.5', '79.5', '7.5', '6.5', '7.5'}, 'ut', {'56.5', '57.5', '58.5', '59.5', ...
    '62.5', '63.5', '64.5', '65.5', '100.5', '101.5', '102.5', '103.5', '13.5', '14.5', '1.5', ...
    '8.5', '12.5', '2.5', '13.5', '80.5', '8.5', '7.5', '8.5'});
% Defining my VOIs after fslmerge not for the exclude_GM VOIs though
%     sn_VOIs = struct('name', {'BStem_voi', 'L_amyg_voi', 'L_cuneus_voi',
%     'L_hipp_voi', 'L_apahc_voi', 'L_ppahc_voi', 'L_PreCG_voi', 'PCC_voi',
%     'PC_voi', ...
%         'R_amyg_voi', 'R_cuneus_voi', 'R_hipp_voi', 'R_apahc_voi',
%         'R_ppahc_voi', 'R_PreCG_voi', 'L_GM_voi', 'L_WM_voi', 'R_GM_voi',
%         'R_WM_voi', 'L_temporal_voi', 'R_temporal_voi'}, ... 'index',
%         {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
%         '12', '13', '14' ,'15', '16', '17', '18', '19'});
sn_VOIs = struct('name', {'BStem_voi', 'L_amyg_voi', 'L_cuneus_voi', 'L_hipp_voi', 'L_apahc_voi', 'L_ppahc_voi', 'L_PreCG_voi', 'R_vcing_voi', 'PCC_voi', 'PC_voi', ...
    'R_amyg_voi', 'R_cuneus_voi', 'R_hipp_voi', 'R_apahc_voi', 'R_ppahc_voi', 'R_PreCG_voi', 'L_vcing_voi', 'L_GM_voi', 'L_WM_voi', 'R_GM_voi', ...
    'R_WM_voi', 'R_temporal_voi', 'L_temporal_voi'}, ...
    'index', {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14' ,'15', '16', '17', '18', '19', '20', '21', '22'});

% These are the VOIs we're using now 

save([dir_main filesep 'sn_VOIs.mat'] , 'sn_VOIs');


%     exclude_VOIs = struct('name', {'L_GM', 'R_GM'}); temporals =
%     struct('name', {'L_temporal', 'R_temporal'});

for i = 1:(size(s_bc_DTI,1))
    DTI_4 = [s_pp_DTI(i).name(1:end-6) 'pp_bc.mif'];
    DTI_4_b0 = [s_bc_DTI(i).name(1:end-9) 'b0.nii.gz'];
    dwi_basename = [s_bc_DTI(i).name(1:end-9)];
    DTI_4_DT = [s_pp_DTI(i).name(1:end-6) 'DT.mif'];
    subj_DTI_out = [DTI_outs(i).folder filesep DTI_outs(i).name];
    DTI_FA_5 = [s_pp_DTI(i).name(1:end-6) 'FA.nii.gz'];
    DTI_ADC_5 = [s_pp_DTI(i).name(1:end-6) 'MD.nii.gz'];
    FA_2_T1 = [s_pp_DTI(i).name(1:end-6) 'FA_2_T1'];
    dh_rf_wm = [s_pp_DTI(i).name(1:end-6) 'Dh_RF_sfwm.txt'];
    to_rf_wm = [s_pp_DTI(i).name(1:end-6) 'To_RF_sfwm.txt'];
    dh_rf_gm = [s_pp_DTI(i).name(1:end-6) 'Dh_RF_gm.txt'];
    dh_rf_csf = [s_pp_DTI(i).name(1:end-6) 'Dh_RF_csf.txt'];
    csd_msmt_wm_fod = [s_pp_DTI(i).name(1:end-6) 'FOD_wm.mif'];
    csd_wm_fod = [s_pp_DTI(i).name(1:end-6) 'csd_FOD_wm.mif'];
    csd_msmt_csf_fod = [s_pp_DTI(i).name(1:end-6) 'FOD_csf.mif'];
    %     upsamp_msmt_wm_fod = [s_pp_DTI(i).name(1:end-6)
    %     'upsamp_csd_FOD_wm.mif']; upsamp_msmt_csf_fod =
    %     [s_pp_DTI(i).name(1:end-6) 'upsamp_csd_FOD_csf.mif'];
    msmt_wm_density = [s_pp_DTI(i).name(1:end-6) 'msmt_wm_density.nii.gz'];
    %     csd_wm_density = [s_pp_DTI(i).name(1:end-6)
    %     'csd_wm_density.nii.gz'];
    dwi_tissues =  [s_pp_DTI(i).name(1:end-6) 'msmt_tissues.nii.gz'];
    wm_fod_corr_mtn  = [s_pp_DTI(i).name(1:end-6) 'mtnc_msmt_wm_fod.mif'];
    csf_fod_corr_mtn  = [s_pp_DTI(i).name(1:end-6) 'mtnc_msmt_csf_fod.mif'];
    %     msmt_wm_peaks_mif = [s_pp_DTI(i).name(1:end-6)
    %     'msmt_wm_peaks.mif'];
    msmt_wm_peaks_nifti = [s_pp_DTI(i).name(1:end-6) 'msmt_wm_peaks.nii.gz'];
    VOIs_4D_inFA = [subj_DTI_out filesep dwi_basename 'FA_VOIs_4D.nii.gz'];
    aal_in_FA = [subj_DTI_out filesep FA_2_T1 '_FA_aal_labels.nii.gz'];
    HOCA_in_FA = [subj_DTI_out filesep FA_2_T1 '_FA_HOCA.nii.gz'];
    HOSCA_in_FA = [subj_DTI_out filesep FA_2_T1 '_FA_HOSCA.nii.gz'];
    MNI_in_FA = [subj_DTI_out filesep FA_2_T1 '_MNI_rl1t.nii.gz'];
    folder = [subj_DTI_out filesep dwi_basename 'VOIs'];
    VOIs_4D = [folder filesep dwi_basename 'VOIs_4D.nii.gz'];
    mkdir(folder);
    
    
    
    % R 04/01/2019 Not resizing anything
    
    % Extracting B0s
    
    unix(['source ~/.bash_profile ; dwiextract -force -bzero -nthreads 7 ' s_bc_DTI(i).folder filesep s_bc_DTI(i).name ' ' s_bc_DTI(i).folder filesep  DTI_4_b0]);
    
    % Creating brain mask of B0s
    
    unix(['source ~/.bash_profile ; antsBrainExtraction.sh -d 3 -a ' s_bc_DTI(i).folder filesep DTI_4_b0 ' -e ' T2_template ' -m ' T2_temp_brain_mask ...
        ' -o ' s_bc_DTI(i).folder filesep dwi_basename 'b0_ -s nii.gz -u 1']);
    
    % dwi2mask is useful for the tensor etc... we can mask later, this is
    % to try and avoid having a rim of artefactually increased FA
    
    unix(['source ~/.bash_profile ; dwi2mask -force -nthreads 7 ' s_bc_DTI(i).folder filesep s_bc_DTI(i).name ' '  s_bc_DTI(i).folder filesep (s_bc_DTI(i).name(1:end-9)) 'mask.nii.gz ']);
    
    % Calculate tensor using ANTs BET mask then generate FA, and MD maps
    
    unix(['source ~/.bash_profile ; dwi2tensor -force -nthreads 7 -mask ' s_bc_DTI(i).folder filesep dwi_basename 'mask.nii.gz ' s_bc_DTI(i).folder filesep DTI_4 ...
        ' ' s_bc_DTI(i).folder filesep DTI_4_DT]);
    
    % output FA and MD images to DTI_output dir
    
    unix(['source ~/.bash_profile ; tensor2metric -force -nthreads 7 -mask ' s_bc_DTI(i).folder filesep dwi_basename 'mask.nii.gz -fa ' subj_DTI_out filesep DTI_FA_5 ' -adc ' subj_DTI_out filesep DTI_ADC_5 ...
        ' ' s_bc_DTI(i).folder filesep DTI_4_DT]);
    
    % calculate RF with Dh
    
    unix(['source ~/.bash_profile ; dwi2response dhollander -force -nthreads 7 -mask ' s_bc_DTI(i).folder filesep dwi_basename 'mask.nii.gz ' s_bc_DTI(i).folder filesep DTI_4 ' ' ...
        s_bc_DTI(i).folder filesep dh_rf_wm ' ' s_bc_DTI(i).folder filesep dh_rf_gm ' ' s_bc_DTI(i).folder filesep dh_rf_csf]);
    
    % Calculate FOD for WM only !?
    
    unix(['source ~/.bash_profile ; dwi2fod -force -nthreads 7 -mask ' s_bc_DTI(i).folder filesep dwi_basename 'mask.nii.gz msmt_csd '  s_bc_DTI(i).folder filesep DTI_4  ' ' ...
        s_bc_DTI(i).folder filesep dh_rf_wm  ' ' s_bc_DTI(i).folder filesep csd_msmt_wm_fod ' ' s_bc_DTI(i).folder filesep dh_rf_csf ' ' s_bc_DTI(i).folder filesep csd_msmt_csf_fod ]);
    
    % should probably do mtnormalise here
    
    unix(['source ~/.bash_profile ; mtnormalise -force -nthreads 7 -mask ' s_bc_DTI(i).folder filesep dwi_basename 'mask.nii.gz  ' ...
        s_bc_DTI(i).folder filesep csd_msmt_wm_fod ' ' s_bc_DTI(i).folder filesep wm_fod_corr_mtn ' ' s_bc_DTI(i).folder filesep csd_msmt_csf_fod ' ' s_bc_DTI(i).folder filesep csf_fod_corr_mtn]);
    
    % upsample FODs to 1.3 isotropic then calculate peaks
    %     unix(['source ~/.bash_profile ; mrresize -force -nthreads 7
    %     -voxel 2 ' s_bc_DTI(i).folder filesep csd_msmt_wm_fod ' '
    %     s_bc_DTI(i).folder filesep upsamp_msmt_wm_fod]); unix(['source
    %     ~/.bash_profile ; mrresize -force -nthreads 7 -voxel 2 '
    %     s_bc_DTI(i).folder filesep csd_msmt_csf_fod ' '
    %     s_bc_DTI(i).folder filesep upsamp_msmt_csf_fod]);
    
    unix(['source ~/.bash_profile ; mrconvert ' s_bc_DTI(i).folder filesep wm_fod_corr_mtn ' ' s_bc_DTI(i).folder filesep msmt_wm_density ' -coord 3 0 ']);
    %     unix(['source ~/.bash_profile ; mrconvert ' s_bc_DTI(i).folder
    %     filesep csd_wm_fod ' ' s_bc_DTI(i).folder filesep csd_wm_density
    %     ' -coord 3 0 ']);
    
    unix(['source ~/.bash_profile ; mrcat ' s_bc_DTI(i).folder filesep msmt_wm_density ' ' s_bc_DTI(i).folder filesep csf_fod_corr_mtn ' ' subj_DTI_out filesep dwi_tissues ' -axis 3']);
    
    % Convert WM FOD to peaks
    unix(['source ~/.bash_profile ; sh2peaks -force -nthreads 7 ' s_bc_DTI(i).folder filesep wm_fod_corr_mtn ' ' s_bc_DTI(i).folder filesep msmt_wm_peaks_nifti]);
    
    % also get a .nii.gz and .bval and .bvec of the peaks for TractSeg
    unix(['source ~/.bash_profile ; mrinfo -force -nthreads 7 ' s_bc_DTI(i).folder filesep s_bc_DTI(i).name ' -export_grad_fsl ' s_bc_DTI(i).folder filesep (msmt_wm_peaks_nifti(1:end-7)) '.bvec ' ...
        s_bc_DTI(i).folder filesep (msmt_wm_peaks_nifti(1:end-7)) '.bval -export_grad_mrtrix ' s_bc_DTI(i).folder filesep (msmt_wm_peaks_nifti(1:end-7)) '_bmatrix.txt'  ]);
    
    
    % Warping MNI 1mm to FA space
    unix(['source ~/.bash_profile ; antsRegistrationSyNQuick.sh -d 3 -n 7 -j 1 -t s -m ' subj_DTI_out filesep DTI_FA_5 ...
        ' -f '  aal_template_brain ' -o ' subj_DTI_out filesep FA_2_T1 'FA_2_MNI  -x ' ...
        subj_DTI_out filesep FA_2_T1 'b0_BrainExtractionMask.nii.gz']);
    
    
    
    % clear previous parallel pool if still present and start a new one
    % with 4 workers
    poolobj = gcp('nocreate');
    delete(poolobj)
    parpool(4);
    
    parfor n = 1:(size(VOIs,2))
        VOI_name = ([folder filesep dwi_basename (VOIs(n).name(7:end))]);
        VOI_lt = (VOIs(n).lt);
        VOI_ut = (VOIs(n).ut);
        
        if (n <= 14)
            labels_file = aal_labels
        elseif (n > 14 & n < 20)
            labels_file = HOSCA
        elseif (n > 19 & n <= 21)
            labels_file = MNI_lobes_rl
        elseif (n > 21 )
            labels_file = JHU_WM
        end
        
        % isolating each VOIs, smoothing, thr at 0.1 and binarizing.
        if (n <= 15 || n >= 22)
            unix(['source ~/.bash_profile ; fslmaths '  labels_file ' -thr ' VOI_lt ' -uthr ' VOI_ut ' -s 1.5 -thr 0.1 -bin ' VOI_name '_inprep.nii.gz']);
        elseif (n > 15 && n <22)
            unix(['source ~/.bash_profile ; fslmaths '  labels_file ' -thr ' VOI_lt ' -uthr ' VOI_ut ' -bin ' VOI_name '_inprep.nii.gz']);
        end
        
    end
    
    unix(['source ~/.bash_profile ; fslmerge -t ' VOIs_4D ' ' folder filesep '*_inprep.nii.gz']);
    
    % This is the correct way of using WIMT in reverse :)
    unix(['source ~/.bash_profile ; WarpTimeSeriesImageMultiTransform 4 ' VOIs_4D ' ' VOIs_4D_inFA ...
        ' -R ' subj_DTI_out filesep DTI_FA_5 ' -i ' subj_DTI_out filesep FA_2_T1 'FA_2_MNI0GenericAffine.mat '  subj_DTI_out filesep FA_2_T1 'FA_2_MNI1InverseWarp.nii.gz --use-NN']);
    
    unix(['source ~/.bash_profile ; WarpImageMultiTransform 3 ' aal_labels ' ' aal_in_FA ...
        ' -R ' subj_DTI_out filesep DTI_FA_5 ' -i ' subj_DTI_out filesep FA_2_T1 'FA_2_MNI0GenericAffine.mat '  subj_DTI_out filesep FA_2_T1 'FA_2_MNI1InverseWarp.nii.gz --use-NN']);
    
    unix(['source ~/.bash_profile ; WarpImageMultiTransform 3 ' HOCA ' ' HOCA_in_FA ...
        ' -R ' subj_DTI_out filesep DTI_FA_5 ' -i ' subj_DTI_out filesep FA_2_T1 'FA_2_MNI0GenericAffine.mat '  subj_DTI_out filesep FA_2_T1 'FA_2_MNI1InverseWarp.nii.gz --use-NN']);
    
    unix(['source ~/.bash_profile ; WarpImageMultiTransform 3 ' HOSCA ' ' HOSCA_in_FA ...
        ' -R ' subj_DTI_out filesep DTI_FA_5 ' -i ' subj_DTI_out filesep FA_2_T1 'FA_2_MNI0GenericAffine.mat '  subj_DTI_out filesep FA_2_T1 'FA_2_MNI1InverseWarp.nii.gz --use-NN']);
    
    unix(['source ~/.bash_profile ; WarpImageMultiTransform 3 ' MNI_lobes_rl ' ' MNI_in_FA ...
        ' -R ' subj_DTI_out filesep DTI_FA_5 ' -i ' subj_DTI_out filesep FA_2_T1 'FA_2_MNI0GenericAffine.mat '  subj_DTI_out filesep FA_2_T1 'FA_2_MNI1InverseWarp.nii.gz --use-NN']);
    
    % The atlases are just for a sanity check, and maybe for later use if
    % necessary, however the VOIs warped to dMRI space in 4D will be used
    % for the tracking.
    
    
    clear n
    
    parfor n = 1:size(sn_VOIs,2)
        
        unix(['source ~/.bash_profile ; fslroi ' VOIs_4D_inFA ' ' folder filesep dwi_basename sn_VOIs(n).name '.nii.gz 0 -1 0 -1 0 -1 ' sn_VOIs(n).index ' 1']);
        
        delete([folder filesep '*_inprep.nii.gz']);
        
    end
end

clear VOI_name

% So now we have exclusionary rois for white matter and grey matter of each
% hemisphere separately with the include/seeding VOIs (those belonging to
% the medial temporal lobe as holes in the exclude VOIs also.. I hope this
% works yo!

% parfor n =1:2
    VOI_name = ([folder filesep dwi_basename]);
    % editing the grey matter exclude VOI
    unix(['source ~/.bash_profile ; fslmaths '  VOI_name sn_VOIs(9).name ' -add ' VOI_name sn_VOIs(4).name ' -add ' VOI_name sn_VOIs(6).name ... 
        ' -add ' VOI_name sn_VOIs(22).name ' -add ' VOI_name sn_VOIs(17).name ' -binv -mul ' VOI_name sn_VOIs(18).name ' ' VOI_name 'L_exclude_GM_voi.nii.gz']);
    unix(['source ~/.bash_profile ; fslmaths '  VOI_name sn_VOIs(9).name ' -add ' VOI_name sn_VOIs(13).name ' -add ' VOI_name sn_VOIs(15).name ...
        ' -add ' VOI_name sn_VOIs(23).name ' -add ' VOI_name sn_VOIs(8).name ' -binv -mul ' VOI_name sn_VOIs(20).name ' ' VOI_name 'R_exclude_GM_voi.nii.gz']);
    % editing the white matter exclude VOI
    unix(['source ~/.bash_profile ; fslmaths '  VOI_name sn_VOIs(9).name ' -add ' VOI_name sn_VOIs(4).name ' -add ' VOI_name sn_VOIs(6).name ... 
        ' -add ' VOI_name sn_VOIs(22).name ' -add ' VOI_name sn_VOIs(17).name ' -add ' VOI_name sn_VOIs(10).name ' -binv -mul ' VOI_name sn_VOIs(19).name ' ' VOI_name 'L_exclude_WM_voi.nii.gz']);
    unix(['source ~/.bash_profile ; fslmaths '  VOI_name sn_VOIs(9).name ' -add ' VOI_name sn_VOIs(13).name ' -add ' VOI_name sn_VOIs(15).name ...
        ' -add ' VOI_name sn_VOIs(23).name ' -add ' VOI_name sn_VOIs(8).name ' -add ' VOI_name sn_VOIs(10).name ' -binv -mul ' VOI_name sn_VOIs(21).name ' ' VOI_name 'R_exclude_WM_voi.nii.gz']);
    % unix(['source ~/.bash_profile ; fslmaths '  VOI_name
    % (sn_VOIs(7).name(1:end)) ' -add ' VOI_name (sn_VOIs(8).name(1:end)) '
    % -add ' VOI_name (sn_VOIs(3).name(1:end)) ...
    %     ' -bin ' VOI_name 'L_Par_rois.nii.gz']);
    % unix(['source ~/.bash_profile ; fslmaths '  VOI_name
    % (sn_VOIs(9).name(1:end)) ' -add ' VOI_name (sn_VOIs(11).name(1:end))
    % ' -add ' VOI_name (sn_VOIs(12).name(1:end)) ...
    %     ' -add ' VOI_name (sn_VOIs(13).name(1:end)) ' -bin ' VOI_name
    %     'R_Temp_rois.nii.gz']);
    % unix(['source ~/.bash_profile ; fslmaths '  VOI_name
    % (sn_VOIs(2).name(1:end)) ' -add ' VOI_name (sn_VOIs(4).name(1:end)) '
    % -add ' VOI_name (sn_VOIs(5).name(1:end)) ...
    %     ' -add ' VOI_name (sn_VOIs(6).name(1:end)) ' -bin ' VOI_name
    %     'L_Temp_rois.nii.gz']);
% end

poolobj = gcp('nocreate');
delete(poolobj)

%% Part 5 ANTs Cortical thickness & AAL labels to native DWI space
s_FA = (dir([dir_DTI_out filesep '*' filesep '*DTI_FA.nii.gz']));
s_DWI = (dir([dir_DTI_int_out filesep 'PT_*_DTI_int_out']));
clear s_T1s;
s_T1s = (dir([dir_t1_int_out filesep '*' filesep '*_T1_bc.nii.gz']));
s_T1_brain = (dir([dir_t1_int_out filesep '*' filesep '*_BrainExtractionBrain.nii.gz']));

% Bringing AAL labels and PD25 template to native T1 space and running ANTs
% Cortical thickness
for i  = 1:(size(s_T1s,1))
    subj_basename = (s_T1s(i).name(1:end-9)); % to create a base name for ANTs CTH output
    dwi_basename = (s_FA(i).name(1:end-9));
    %     Cort_th_subj = [dir_t1_out filesep subj_basename 'output' filesep
    %     subj_basename 'cortical_thickness' filesep subj_basename
    %     'antsCTh_'];
    Cort_th_subj = [dir_t1_out filesep subj_basename 'output' filesep  subj_basename 'cortical_thickness'];
    mkdir(Cort_th_subj);
    tpms = (dir([dir_t1_out filesep subj_basename 'output' filesep  subj_basename 'cortical_thickness' filesep subj_basename 'BrainSegmentationPosteriors*.nii.gz']));
    empty_tpm = [dir_t1_out filesep subj_basename 'output' filesep  subj_basename 'cortical_thickness' filesep subj_basename 'empty_tpm.nii.gz'];
    p5tt_4d = [dir_t1_out filesep subj_basename 'output' filesep  subj_basename 'cortical_thickness' filesep subj_basename 'p5tt_4D.nii.gz'];
    p5tt_4d_in_FA = [s_FA(i).folder filesep  dwi_basename 'p5tt_4D_inFA.nii.gz'];
    aal_in_subj = [s_T1s(i).folder filesep 'aal_labels_in_' subj_basename 'native.nii.gz'];
    sub_DTI_int_out = [s_DWI(i).folder filesep s_DWI(i).name];
    
    unix(['source ~/.bash_profile ; antsRegistrationSyNQuick.sh -d 3 -n 7 -j 1 -t s -m ' T1_template_brain ...
        ' -f ' s_T1s(i).folder filesep subj_basename 'BrainExtractionBrain.nii.gz -o ' ...
        s_T1s(i).folder filesep 'PD25brain_2_' subj_basename 'brain  -x ' ...
        s_T1s(i).folder filesep subj_basename 'BrainExtractionMask.nii.gz']);
    
    unix(['source ~/.bash_profile ; antsRegistrationSyNQuick.sh -d 3 -n 7 -j 1 -t s -m ' aal_template_brain ...
        ' -f ' s_T1s(i).folder filesep subj_basename 'BrainExtractionBrain.nii.gz -o ' ...
        s_T1s(i).folder filesep 'MNI1mbrain_2_' subj_basename 'brain  -x ' ...
        s_T1s(i).folder filesep subj_basename 'BrainExtractionMask.nii.gz']);
    
    unix(['source ~/.bash_profile ; WarpImageMultiTransform 3 ' aal_labels ' ' aal_in_subj ...
        ' -R ' s_T1s(i).folder filesep subj_basename 'BrainExtractionBrain.nii.gz ' ...
        s_T1s(i).folder  filesep 'MNI1mbrain_2_' subj_basename 'brain1Warp.nii.gz ' s_T1s(i).folder filesep 'MNI1mbrain_2_' subj_basename 'brain0GenericAffine.mat --use-NN']); % still using NN as ML is rubbish
    
    % extracting cortical thickness maps
    unix(['source ~/.bash_profile ; antsCorticalThickness.sh -d 3 -a ' s_T1s(i).folder  filesep s_T1s(i).name  ' -e ' T1_template ' -m ' T1_temp_brain_mask ...
        ' -p ' dir_templates filesep 'antsCTh_PD25_priors_0%d.nii.gz ' ' -r '  aal_in_subj ' -k 1 -g 1 -q 1 -o ' Cort_th_subj filesep subj_basename]);
    
    % Prepping the T1 tpms
    unix(['source ~/.bash_profile ; fslmaths ' tpms(1).folder filesep tpms(1).name ' -mul 0 ' empty_tpm ]);
    
    unix(['source ~/.bash_profile ; fslmerge -t ' p5tt_4d ' ' tpms(1).folder filesep tpms(2).name ' ' tpms(1).folder filesep tpms(4).name ...
        ' ' tpms(1).folder filesep tpms(3).name ' ' tpms(1).folder filesep tpms(1).name ' ' empty_tpm]);
    
    unix(['source ~/.bash_profile ; antsRegistrationSyNQuick.sh -d 3 -n 7 -j 1 -t s -m ' s_FA(i).folder filesep s_FA(i).name ...
        ' -f ' s_T1_brain(i).folder filesep s_T1_brain(i).name ' -o ' sub_DTI_int_out filesep dwi_basename 'FA_2_T1_brain -x ' sub_DTI_int_out filesep dwi_basename 'b0_BrainExtractionMask.nii.gz']);
    
    % To get your tpms into DWI space - this finally works now
    unix(['source ~/.bash_profile ; WarpTimeSeriesImageMultiTransform 4 ' p5tt_4d ' ' p5tt_4d_in_FA ' -R ' s_FA(i).folder filesep s_FA(i).name ' -i ' ...
        sub_DTI_int_out filesep dwi_basename 'FA_2_T1_brain0GenericAffine.mat ' sub_DTI_int_out filesep dwi_basename 'FA_2_T1_brain1InverseWarp.nii.gz']);
    
    % Then just sanity check with 5ttcheck
    unix(['source ~/.bash_profile ; 5ttcheck -masks ' sub_DTI_int_out filesep dwi_basename '5tt_fail.nii.gz ' p5tt_4d_in_FA]);
    % Consider also adding FreeSurfer
end


%% Part 6 Fiber tracking and tckmap
% First we need to define a few variables for the tckgen command Why not
% use ACT here, just need to warp, that is all.
save ([dir_main filesep 'ws_pre_ft.mat']);
clear all;
clc;

dir_main   = '/Volumes/LaCie/MTLE_HS';
dir_source_DTI = [dir_main filesep 'NII' filesep 'DTI' filesep 'orig'];
dir_source_T1 = [dir_main filesep 'NII' filesep 'T1s' filesep 'orig'];
dir_t1_int_out = [dir_main filesep 'NII' filesep 'T1s' filesep 'intermediate'];
dir_t1_out = [dir_main filesep 'NII' filesep 'T1s' filesep 'T1_output'];
% dir_t1_temp_wskull = [dir_t1_int_out filesep 'T1_template_wskull'];
% dir_t1_temp_brain = [dir_t1_int_out filesep 'T1_template_brain'];
dir_DTI_int_out = [dir_main filesep 'NII' filesep 'DTI' filesep 'intermediate'];
dir_DTI_out = [dir_main filesep 'NII' filesep 'DTI' filesep 'DTI_output'];
dir_log    = [dir_main filesep 'LOG'];
dir_templates = [dir_main filesep 'templates'];
tck_jobs = dir([dir_DTI_out filesep '*']);

% need to define a nice variable containing all my VOIs.... then try
% tracking again :) R - 05/01/2018
% also need to add in MPFC for dorsal cingulum.

% we need the peaks, the FAs, the fods, the masks, the dwis, and the p5tt  + the VOIs
s_peaks = (dir([dir_DTI_int_out filesep '*' filesep '*peaks.nii.gz']));
FAs = (dir([dir_DTI_out filesep '*' filesep '*_DTI_FA.nii.gz']));
dwi_p5tt = (dir([dir_DTI_out filesep '*' filesep '*p5tt_4D_inFA.nii.gz']));
s_fods = (dir([dir_DTI_int_out filesep '*' filesep '*mtnc_msmt_wm_fod.mif']));
s_dwis = (dir([dir_DTI_int_out filesep '*' filesep '*_pp_bc.mif']));
s_masks = (dir([dir_DTI_int_out filesep '*' filesep '*BrainExtractionMask.nii.gz']));

VOI_dirs = struct(dir([dir_DTI_out filesep '*' filesep '*VOIs']));
VOIs_comp = struct([]);

parfor i = 1:size(VOI_dirs,1);
  
    
    VOIs_comp(i).name = (VOI_dirs(i).name(1:end-9));
    VOIs = (dir([VOI_dirs(i).folder filesep VOI_dirs(i).name filesep '*voi.nii.gz']));
    VOIs_comp(i).voi_name = ({VOIs.name}');
    VOIs_comp(i).voi_folder = ([VOI_dirs(i).folder filesep VOI_dirs(i).name]);

end


% Loop over the peaks files and do tckgen with include and exclude rois
% defined above, if using ACT need to implement warping of tpms from t1 to
% b0 space... should be simple enought, perhaps another 80 lines.


% This isn't working yet... found a stupid bug from yesterday, should try
% it again with FACT now that this bug is out. result from Tensor_Prob
% looks quite good (for whole brain) with 1 mil

for i = 1:(size(s_peaks,1));
    basename = (s_peaks(i).name(1:end-20));
    int_output = s_peaks(i).folder;
    final_output = FAs(i).folder;
    tck_out = ([final_output filesep 'tck_output']);
    mkdir(tck_out);
    
    % whole brain then tckedit
    unix(['source ~/.bash_profile ; tckgen -algorithm Tensor_Det -force -nthreads 14 -select 5000000 ' int_output filesep s_dwis(i).name ' ' tck_out filesep basename 'TD_wholebrain.tck  -seed_dynamic ' ...
        int_output filesep s_fods(i).name   ' -mask ' int_output filesep s_masks(i).name  ' -act ' final_output filesep dwi_p5tt(i).name]);
%         ' -include ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{13} ...
%         ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{1} ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{8} ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{17}]);
    
    
    unix(['source ~/.bash_profile ; tckedit -force -nthreads 14 -select 10000 ' tck_out filesep basename 'TD_wholebrain.tck ' tck_out filesep basename 'TP_L_phC.tck ' ...
        ' -mask ' int_output filesep s_masks(i).name ...
%         ' -include ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{13} ...
        ' -include ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{14} ' -include ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{15} ...
        ' -include ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{5} ...
        ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{1} ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{8} ...
        ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{9} ...
        ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{18} ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{16}]);
    
    unix(['source ~/.bash_profile ; tckgen -algorithm FACT -force -nthreads 7 -select 5000 ' int_output filesep s_peaks(i).name ' ' final_output filesep basename 'FACT_R_phC.tck  -seed_gmwmi ' ...
        VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{18} ' -mask ' int_output filesep s_masks(i).name  ' -act ' final_output filesep dwi_p5tt(i).name ...
        ' -include ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{13} ...
        ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{1} ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{21} ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{4}]);
    
    
    
    unix(['source ~/.bash_profile ; tckgen -algorithm iFOD2 -force -nthreads 7 -select 5000 ' int_output filesep s_fods(i).name ' ' final_output filesep basename 'TP_L_CST.tck  -seed_image ' ...
        VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{1}  ' -act ' final_output filesep dwi_p5tt(i).name ' -mask ' int_output filesep s_masks(i).name  ...
        ' -include ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{3} ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{17}]);
    
    % ' -include ' VOIs_comp(i).voi_folder filesep
    % VOIs_comp(i).voi_name{13} ... ' -exclude ' VOIs_comp(i).voi_folder
    % ' -act ' final_output filesep dwi_p5tt(i).name ' -include ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{12} ' ' filesep VOIs_comp(i).voi_name{1} ' -exclude ' VOIs_comp(i).voi_folder filesep VOIs_comp(i).voi_name{8} 
    
    
    
    
    
    
    
    %     VOIs_dir = s_Bstem(i).folder;
    
    % consider using a -maxlength of 200 mm, this may be set to even lower
    % if looking for short connections specifically and also using
    % -minlength may help clean up the tracts a bit, crude tracking first
    % %followed by some stats on to derive reasonable case specific min and
    % max length or simple a hard coded one, as an optional input...
    
    
    % The left phCing with the 4 methods of tracking
    %     unix(['source ~/.bash_profile ; tckgen -algorithm Tensor_Prob
    %     -force -nthreads 7 -select 10000 ' int_output filesep
    %     s_dwis(i).name ' ' final_output filesep basename 'TP_L_phC_TP.tck
    %     -seed_dynamic ' ...
    %         int_output filesep s_fods(i).name  ' -act ' final_output
    %         filesep dwi_p5tt(i).name ' -mask ' int_output filesep
    %         s_masks(i).name ' -include ' VOIs_dir filesep
    %         s_l_hipp(i).name ... ' -include ' VOIs_dir filesep
    %         s_L_par_rois(i).name ' -exclude ' VOIs_dir filesep
    %         s_R_wm(i).name ' -exclude ' VOIs_dir filesep s_Bstem(i).name
    %         ]);
    
    % need to try this with -seed_image and use the hippocampus as the
    % seed.
    
    % perhaps I need to rethink this... why not track the whole cingulum
    % bundle!! use the MPFC, PCC and hippos + apaHC
    unix(['source ~/.bash_profile ; tckgen -algorithm Tensor_Prob -force -nthreads 7 -select 5000 ' int_output filesep s_dwis(i).name ' ' final_output filesep basename 'TP_L_phC_TP.tck  -seed_gmwmi ' ...
        VOIs_dir filesep s_l_hipp(i).name ' -act ' final_output filesep dwi_p5tt(i).name ' -mask ' int_output filesep s_masks(i).name  ...
        ' -include ' VOIs_dir filesep s_L_par_rois(i).name ' -exclude ' VOIs_dir filesep s_R_wm(i).name ' -exclude ' VOIs_dir filesep s_Bstem(i).name ]);
    % and once for the right side
    unix(['source ~/.bash_profile ; tckgen -algorithm Tensor_Prob -force -nthreads 7 -select 5000 ' int_output filesep s_dwis(i).name ' ' final_output filesep basename 'TP_R_phC_TP.tck  -seed_gmwmi ' ...
        VOIs_dir filesep s_r_hipp(i).name ' -act ' final_output filesep dwi_p5tt(i).name ' -mask ' int_output filesep s_masks(i).name  ...
        ' -include ' VOIs_dir filesep s_R_par_rois(i).name ' -exclude ' VOIs_dir filesep s_L_wm(i).name ' -exclude ' VOIs_dir filesep s_Bstem(i).name ]);
    
    unix(['source ~/.bash_profile ; tckgen -algorithm Tensor_Det -force -nthreads 7 -select 10000 ' int_output filesep s_dwis(i).name ' ' final_output filesep basename 'TP_L_phC_TD.tck  -seed_dynamic ' ...
        int_output filesep s_fods(i).name  ' -act ' final_output filesep dwi_p5tt(i).name ' -mask ' int_output filesep s_masks(i).name ' -include ' VOIs_dir filesep s_l_hipp(i).name ...
        ' -include ' VOIs_dir filesep s_L_par_rois(i).name ' -exclude ' VOIs_dir filesep s_R_wm(i).name ' -exclude ' VOIs_dir filesep s_Bstem(i).name ]);
    unix(['source ~/.bash_profile ; tckgen -algorithm FACT -force -nthreads 7 -select 10000 ' int_output filesep s_peaks(i).name ' ' final_output filesep basename 'TP_L_phC_FACT.tck  -seed_dynamic ' ...
        int_output filesep s_fods(i).name  ' -act ' final_output filesep dwi_p5tt(i).name ' -mask ' int_output filesep s_masks(i).name ' -include ' VOIs_dir filesep s_l_hipp(i).name ...
        ' -include ' VOIs_dir filesep s_L_par_rois(i).name ' -exclude ' VOIs_dir filesep s_R_wm(i).name ' -exclude ' VOIs_dir filesep s_Bstem(i).name ]);
    unix(['source ~/.bash_profile ; tckgen -algorithm iFOD2 -force -nthreads 7 -select 10000 ' int_output filesep s_fods(i).name ' ' final_output filesep basename 'TP_L_phC_iFOD2.tck  -seed_dynamic ' ...
        int_output filesep s_fods(i).name  ' -act ' final_output filesep dwi_p5tt(i).name ' -mask ' int_output filesep s_masks(i).name ' -include ' VOIs_dir filesep s_l_hipp(i).name ...
        ' -include ' VOIs_dir filesep s_L_par_rois(i).name ' -exclude ' VOIs_dir filesep s_R_wm(i).name ' -exclude ' VOIs_dir filesep s_Bstem(i).name ]);
    
    % The Right phCing with the 4 methods of tracking
    unix(['source ~/.bash_profile ; tckgen -algorithm Tensor_Prob -force -nthreads 7 -select 10000 ' int_output filesep s_dwis(i).name ' ' final_output filesep basename 'TP_R_phC_TP.tck  -seed_dynamic ' ...
        int_output filesep s_fods(i).name  ' -act ' final_output filesep dwi_p5tt(i).name ' -mask ' int_output filesep s_masks(i).name ' -include ' VOIs_dir filesep s_r_hipp(i).name ...
        ' -include ' VOIs_dir filesep s_R_par_rois(i).name ' -exclude ' VOIs_dir filesep s_L_wm(i).name ' -exclude ' VOIs_dir filesep s_Bstem(i).name ])
    unix(['source ~/.bash_profile ; tckgen -algorithm Tensor_Det -force -nthreads 7 -select 10000 ' int_output filesep s_dwis(i).name ' ' final_output filesep basename 'TP_R_phC_TD.tck  -seed_dynamic ' ...
        int_output filesep s_fods(i).name  ' -act ' final_output filesep dwi_p5tt(i).name ' -mask ' int_output filesep s_masks(i).name ' -include ' VOIs_dir filesep s_r_hipp(i).name ...
        ' -include ' VOIs_dir filesep s_R_par_rois(i).name ' -exclude ' VOIs_dir filesep s_L_wm(i).name ' -exclude ' VOIs_dir filesep s_Bstem(i).name ]);
    unix(['source ~/.bash_profile ; tckgen -algorithm FACT -force -nthreads 7 -select 10000 ' int_output filesep s_peaks(i).name ' ' final_output filesep basename 'TP_R_phC_FACT.tck  -seed_dynamic ' ...
        int_output filesep s_fods(i).name  ' -act ' final_output filesep dwi_p5tt(i).name ' -mask ' int_output filesep s_masks(i).name ' -include ' VOIs_dir filesep s_r_hipp(i).name ...
        ' -include ' VOIs_dir filesep s_R_par_rois(i).name ' -exclude ' VOIs_dir filesep s_L_wm(i).name ' -exclude ' VOIs_dir filesep s_Bstem(i).name ]);
    unix(['source ~/.bash_profile ; tckgen -algorithm iFOD2 -force -nthreads 7 -select 10000 ' int_output filesep s_fods(i).name ' ' final_output filesep basename 'TP_R_phC_iFOD2.tck  -seed_dynamic ' ...
        int_output filesep s_fods(i).name  ' -act ' final_output filesep dwi_p5tt(i).name ' -mask ' int_output filesep s_masks(i).name ' -include ' VOIs_dir filesep s_r_hipp(i).name ...
        ' -include ' VOIs_dir filesep s_R_par_rois(i).name ' -exclude ' VOIs_dir filesep s_L_wm(i).name ' -exclude ' VOIs_dir filesep s_Bstem(i).name ]);
    
    
    %     unix(['source ~/.bash_profile ; tckmap -force -nthreads 7 ']);
end

%TractSeg on peaks!!!!

toc
diary off
