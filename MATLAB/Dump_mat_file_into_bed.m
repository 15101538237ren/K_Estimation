function Dump_mat_file_into_bed()
%Load data
MAT_FILE_PATH = 'DATA/AllDatReps123_Chr1WT.mat';
load(MAT_FILE_PATH, 'AllDat', 'sites');
%offset the cpg site by +1
sites = sites + 1;


% Processing Methylated data
OUT_METHY_FILE_PATH = 'DATA/AllDatReps123_Chr1WT_METHY.bed';

METHY_DATA = zeros(length(sites), 5);
METHY_DATA(:, 1) = sites;
METHY_DATA(:, 2 : end)  = AllDat(:,[1, 2, 4, 6], 1);

fileID = fopen(OUT_METHY_FILE_PATH, 'w');
fprintf(fileID,'%d\t%d\t%d\t%d\t%d\n',METHY_DATA');
fclose(fileID);

% Processing Unmethylated data
OUT_UNMETHY_FILE_PATH = 'DATA/AllDatReps123_Chr1WT_UNMETHY.bed';
UNMETHY_DATA = zeros(length(sites), 5);
UNMETHY_DATA(:, 1) = sites;
UNMETHY_DATA(:, 2 : end)  = AllDat(:,[1, 2, 4, 6], 2);

fileID = fopen(OUT_UNMETHY_FILE_PATH, 'w');
fprintf(fileID,'%d\t%d\t%d\t%d\t%d\n',UNMETHY_DATA');
fclose(fileID);
end