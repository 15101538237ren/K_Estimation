clear;
close all;
clc;
comb_size = 1;
chr_size = 1;
DATA_DIR = '../DATA/Repli_BS/MATLAB_DATA/';

OUT_DIR = '../DATA/Repli_BS/K_RATES/';
if ~exist(OUT_DIR)
    mkdir(OUT_DIR);
end
comb_target = [1, 41];
for i = 1 :length(comb_target)
    for j = 1 : chr_size
       comb_i = num2str(comb_target(i));
       data_fp = strcat(DATA_DIR, comb_i, '/', 'chr', num2str(j) , '.mat');
       
       OUT_SUB_DIR = strcat(OUT_DIR, comb_i);
       if ~exist(OUT_SUB_DIR)
            mkdir(OUT_SUB_DIR);
       end
       out_fp = strcat(OUT_DIR, comb_i, '/', 'chr', num2str(j) , '.mat');
       FitMLRates_Protocol1a_RandomT(data_fp, out_fp);
    end
end