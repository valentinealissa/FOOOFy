                            %% DBS EPHYS DATA QC %%
% Last update: 4/25/22
% Alissa Valentine

                            %%% CONTENT %%%                          
% This script produces a .csv/.mat with files to be used for processing based on
% recording QC measures 

% I. Input files
% II. Output files
% III. QC files
% IV. Save

close all
clear

%% 1. Input files
indir = '/Users/alissavalentine/ephys_rotation/';

% changing to input file directory
cd(indir);

% loading CSV
infile = load('depth_duration.mat').mat_file;

%% 2. Output files
outfile = 'qc_depth_duration.mat';

%% 3. QC files

% creating empty array to hold future values
final_array = ["patient_id", "region_side", "depth_mm", "duration_sec", "file_path"];

% QC electrode depth and duration
n = length(infile);
for i = 1:n
    
    if startsWith(infile(i,2), "STN")
        % STN depth range: -5 to 0
        depth = str2double(infile(i,3));
        % recording duration: >5 sec
        dur = str2double(infile(i,4));
        if (dur >= 5)
        %if (-5 < depth) && (depth < 0) && (dur >= 5)
            final_array = cat(1, final_array, infile(i,:));
        end
    end

    if startsWith(infile(i,2), "VIM")
        % VIM depth range: -4 to 0
        depth = str2double(infile(i,3));
        % recording duration: >5 sec
        dur = str2double(infile(i,4));
        if (dur >= 5)
        %if (-4 < depth) && (depth < 0) && (dur >= 5)
            final_array = cat(1, final_array, infile(i,:));
        end
    end
    
    if startsWith(infile(i,2), "GPi")
        % GPi depth range: -6 to 0
        depth = str2double(infile(i,3));
        % recording duration: >5 sec
        dur = str2double(infile(i,4));
        if (dur >= 5)
        %if (-6 < depth) && (depth < 0) && (dur >= 5)
            final_array = cat(1, final_array, infile(i,:));
        end
    end

end

%% 4. Save
save(outfile, 'final_array');
writematrix(final_array,'qc_depth_duration.csv');
