function Plot_correlation_for_all_chromosomes()
    chr_size = 22;
    K_1_BASE_PATH = '../DATA/Repli_BS/K_RATES/1/';
    K_41_BASE_PATH = '../DATA/Repli_BS/K_RATES/41/';
   OUT_FIG_DIR = 'Figures/EXPERIMENT_1_and_41/';
   if ~exist(OUT_FIG_DIR)
        mkdir(OUT_FIG_DIR);
   end
   for i = 1 : chr_size
       K1_PATH = strcat(K_1_BASE_PATH, 'chr', num2str(i) , '.mat');
       K41_PATH = strcat(K_41_BASE_PATH, 'chr', num2str(i) , '.mat');
       Figure_path = strcat(OUT_FIG_DIR, 'chr', num2str(i) , '.pdf');
       Plot_Correlation(K1_PATH, K41_PATH, Figure_path);
   end
end