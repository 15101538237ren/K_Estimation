function Fit_K_for_all_combination_replicates()
clear;
close all;
clc;
comb_size = 41;
chr_size = 1;
DATA_DIR = '../DATA/Repli_BS/MATLAB_DATA/';

OUT_DIR = '../DATA/Repli_BS/K_RATES_NON_RANDOM/';
if ~exist(OUT_DIR)
    mkdir(OUT_DIR);
end

for i = comb_size : comb_size
    for j = 1 : chr_size
       progress=strcat('combine idx: ', num2str(i), ' chr', num2str(j));
       disp(progress);
       comb_i = num2str(i);
       data_fp = strcat(DATA_DIR, comb_i, '/', 'chr', num2str(j) , '.mat');
       
       OUT_SUB_DIR = strcat(OUT_DIR, comb_i);
       if ~exist(OUT_SUB_DIR)
            mkdir(OUT_SUB_DIR);
       end
       out_fp = strcat(OUT_DIR, comb_i, '/', 'chr', num2str(j) , '.mat');
       FitMLRates_Protocol1a(data_fp, out_fp);
    end
end
end
