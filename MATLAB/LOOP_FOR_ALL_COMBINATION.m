comb_size = 1;
chr_size = 1;
DATA_DIR = '../DATA/Repli_BS/MATLAB_DATA/';

OUT_DIR = '../DATA/Repli_BS/K_RATES/';
if ~exist(OUT_DIR)
    mkdir(OUT_DIR);
end

for i= 1 : comb_size
    for j = 1 : chr_size
       data_fp = strcat(DATA_DIR, num2str(i), '/', 'chr', num2str(j) , '.mat');
       
       OUT_SUB_DIR = strcat(OUT_DIR, num2str(i));
       if ~exist(OUT_SUB_DIR)
            mkdir(OUT_SUB_DIR);
       end
       
       out_fp = strcat(OUT_DIR, num2str(i), '/', 'chr', num2str(j) , '.mat');
       
       FitMLRates_Protocol_honglei(data_fp, out_fp)
    end
end