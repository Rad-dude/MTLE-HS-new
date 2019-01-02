% MTLE-HS-dcm2nii
% v1.0 - dd 06/10/2018
% @ radwan - sunaert

% Configure my stuff...
% the script fails for subjects 4 and 13, not sure why (probably related to
% search string.. but eff it
% now also getting rid of spaces in the subfold names
clear all
clc
n_cpu = 4;
dir_main   = '/Users/aradwa0/MR-data/MTLE_HS';


% Other dirs
dir_source = [dir_main filesep 'sorted_dicoms'];
dir_output_DTI = [dir_main filesep 'NII' filesep 'DTI'];
dir_orig_DTI_nii = [dir_output_DTI filesep 'orig'];
dir_output_T1 = [dir_main filesep 'NII' filesep 'T1s'];
dir_T1_orig = [dir_output_T1 filesep 'orig'];
dir_log    = [dir_main filesep 'LOG'];

% Say hello
str = sprintf('\n%s\n\n','*** cp_convert_DICOM_to_NII');
disp(str);

% Create folders
mkdir(dir_log); % create log folder
cmd = ['echo > ' dir_log filesep 'cp_convert_DICOM_to_NII.log'];
unix(cmd);
mkdir(dir_output_DTI); % making output dirs
mkdir(dir_orig_DTI_nii);
mkdir(dir_output_T1);
mkdir(dir_T1_orig);
diary([dir_log filesep 'command_log.txt']); % creating command log


% Search the source for subjects
cd(dir_source);
s = dir;
ns = size(s,1)-3;

% Cycle over subjects & convert each of them

for i = 1:ns
    dir_source_dcm = [dir_source filesep s(i+3).name filesep];
    if i < 10
        pre = '00';
    elseif i < 100;
        pre = '0';
    end
    subj = ['PT_' pre num2str(i)];
    
    disp(['*** NOW @ ' subj]);
    
    %     % Convert the DTI
    subfold = (dir(dir_source_dcm));
    
    strcmpstatus = strcmp(subfold(3).name, '.DS_Store');
    if strcmpstatus == 0,
        subfold = subfold(3).name;
    else
        try
            subfold = subfold(4).name;
        catch
        end
    end
    %     .DS_store and replacing with subfold(4).name.
    A = strrep(subfold,' - ','-');
    unix(['mv -f ' dir_source_dcm '*_* ' dir_source_dcm A])
    subfold = A;
    TFE = (dir([dir_source_dcm filesep subfold filesep '*TFE*']));% find a way to pick out cells containing tfe in name
    DTI = (dir([dir_source_dcm filesep subfold filesep '*DTI*']));% find a way to pick cells containing DTI in name
    str = sprintf('\n%s\n\n','------  Converting the DTIs -----------------------------------------');disp(str);
    %     subj_out = [dir_output_DTI filesep subj '_DTI'];
    subj_out = [subj '_DTI'];
    %     cmd = ['mrconvert -quiet -nthreads ' num2str(n_cpu) ' -force -export_grad_fsl ' subj_out '.bvec '  subj_out '.bval ' dir_source_dcm ' ' subj_out '.nii'];
    %     unix(cmd);
    for d=1:size(DTI,1);
        [status, ~] = unix(['ls ' subj_out '.nii']);
        if status == 1
            %             [mrstatus, ~] = unix(['mrconvert -quiet -nthreads ' num2str(n_cpu) ' -export_grad_fsl ' subj_out '.bvec '  subj_out '.bval ' ...
            %                 ' -export_grad_mrtrix ' subj_out '_bmatrix.txt ' dir_source_dcm subfold filesep DTI(d).name ' ' subj_out '.nii']);
            [~, bla] = unix(['dcm2niix -z y -o ' dir_output_DTI ' -f ' subj_out ' ' dir_source_dcm subfold filesep DTI(d).name])
            k = strfind(bla,'No valid DICOM files were found');
            if k == 163;
                disp(['This folder ' DTI(d).name ' has no useful dicom data, moving to next folder.']);
            else
                disp(['This folder ' DTI(d).name ' has useful dicom data, good, DTI for ' subj ' are converted.']);
            end
            [status, ~] = unix(['ls ' subj_out '.nii']);
            if status == 0;
                break
            end
        end
    end
    
    % Convert the T1
    str = sprintf('\n%s\n\n','------  Converting the T1WIs -----------------------------------------');disp(str);
    %     subj_out = [dir_output_T1 filesep subj '_DTI_T1_orig'];
    %     subj_out2 = [dir_output_T1 filesep subj '_DTI_T1_std'];
    subj_out = [subj '_DTI_T1_orig'];
    for t=1:size(TFE,1);
        [status, ~] = unix(['ls ' subj_out '.nii']);
        if status == 1
            %             [mrstatus, ~] = unix(['mrconvert -quiet -nthreads ' num2str(n_cpu) ' -export_grad_fsl ' subj_out '.bvec '  subj_out '.bval ' ...
            %                 ' -export_grad_mrtrix ' subj_out '_bmatrix.txt ' dir_source_dcm subfold filesep DTI(d).name ' ' subj_out '.nii']);
            [~, bla] = unix(['dcm2niix -z y -o ' dir_output_T1 ' -f ' subj_out ' ' dir_source_dcm subfold filesep TFE(t).name])
            k = strfind(bla,'No valid DICOM files were found');
            if k == 163;
                disp(['This folder ' TFE(t).name ' has no useful dicom data, moving to next folder.']);
            else
                disp(['This folder ' TFE(t).name ' has useful dicom data, good, T1WIs for ' subj ' are converted.']);
            end
            [status, ~] = unix(['ls ' subj_out '.nii']);
            if status == 0;
                break
            end
        end
    end
    %             if using mrconvert, use the lines below..
    %             [mrstatus, ~] = unix(['mrconvert -quiet -nthreads ' num2str(n_cpu) ' ' ...
    %                 dir_source_dcm subfold filesep TFE(t).name ' ' subj_out '.nii']);
    %             if mrstatus == 1
    %                 disp(['This folder ' TFE(t).name ' has no useful dicom data, moving to next folder.']);
    %             else
    %                 disp(['This folder ' TFE(t).name ' has useful dicom data, good, T1WIs for ' subj ' are converted.']);
    %             end
    %             [status, ~] = unix(['ls ' subj_out '.nii']);
    %             if status == 0;
    %                 [reori, ~]= unix(['fslreorient2std ' subj_out '.nii' ' ' subj_out2]);
    %                 unix(['gunzip ' dir_output_T1 filesep '*.gz']);
    %                 break
    %             end
    %     end
    
    cmd = ['echo ''' subj ' - ' s(i+3).name ''''  ' >> ' dir_log filesep 'cp_convert_DICOM_to_NII.log'];
    unix(cmd);
    
end

unix(['mv ' dir_output_DTI filesep '*DTI.* ' dir_orig_DTI_nii]);
unix(['mv ' dir_output_T1 filesep '*DTI_T1_orig* ' dir_T1_orig]);

diary off
