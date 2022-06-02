                            %% DBS EPHYS DATA AGGREGATION %%
% Last update: 4/25/22
% Alissa Valentine

                            %%% CONTENT %%%                          
% This script produces a .csv and .mat file containing the following descriptive stats:
    % patient id
    % region and side
    % electrode depth (mm)
    % recording duration (sec)
    % recording file path

% I. Input files
% II. Output files
% III. Parse files
% IV. Save

close all
clear

%% 1. Input files
indir = '/Users/alissavalentine/ephys_rotation/MATfiles/';
%indir = '/Users/alissavalentine/Documents/MATLAB/dbs_ephys/data/'; 
% update path for wherever your files are! 

%% 2. Output files
outdir = '/Users/alissavalentine/ephys_rotation/';
% update path for wherever you want the CSV file to be saved!

%% 3. Parse JSON
% changing to input file directory
cd(indir);
% collecting relevant directory file info
infiles = dir('PT-*');

% creating empty array to hold future values in .csv
final_array = ["patient_id", "region_side", "depth_mm", "duration_sec", "file_path"];

% creating empty array to hold future values in .mat
mat_file = [];

n = length(infiles);
for i = 1:n
        
    % changing directory to one of the patient's files
    pt_dir = string(indir) + string(infiles(i).name);
    cd(pt_dir);
    
    % counting brain regions
    regions.STNLeft = dir('STNL*');
    regions.STNRight = dir('STNR*');
    regions.VIMLeft = dir('VIML*');
    regions.VIMRight = dir('VIMR*');
    regions.GPiLeft = dir('GPiL*');
    regions.GPiRight = dir('GPiR*');
    region_names = fieldnames(regions);
    
    % saving file info per brain region
    for j = 1:length(region_names)
        r = string(region_names(j));
        
        % creating column for brain region
        if ~isempty(regions.(region_names{j}))
            
            % finding file names
            region_dir = pt_dir + '/' + string(r) + '/';
            cd(region_dir);
            files = dir('*.mat');
            
            % iterating through each file
            for k = 1:length(files)

                % getting rid of the ALL file
                idx = strfind(files(k).name,'ALL');
                if isempty(idx)
                    files(k).name;
                    file_name = strsplit(files(k).name, {'_'});
                    
                    % extracting electrode depth value
                    d = strsplit(string(file_name(6)), {'.mat'});
                    d = extractAfter(d(1),1);
                    
                    % extracting recording duration value
                    data = load(files(k).name);
                    data_name = string(fieldnames(data));
                    srate = data.(data_name).ContinuousChannels.ADFrequency;
                    fragments = data.(data_name).ContinuousChannels.Fragments;
                    time = double(fragments)/double(srate);
                    
                    % add the file path to column
                    path = files(k).folder + "/" + files(k).name;
                    
                    % saving values in stats table
                    stat = [string(infiles(i).name) r d time path];
                    final_array = cat(1, final_array, stat);
                    mat_file = cat(1, mat_file, stat);
                    
                    % saving values in mat file table
                    %stat( = [string(infiles(i).name) r d time path data.(data_name)];
%                     mat_file(i,1) = string(infiles(i).name);
%                     mat_file(i,2) = r;
%                     mat_file(i,3) = d;
%                     mat_file(i,4) = time;
%                     mat_file(i,5) = path;
%                     mat_file(i,6) = data.(data_name);
                end
            end
        end
    end 
end

%% 4. Save
cd(outdir);
save('depth_duration.mat', 'mat_file');
writematrix(final_array,'depth_duration.csv');