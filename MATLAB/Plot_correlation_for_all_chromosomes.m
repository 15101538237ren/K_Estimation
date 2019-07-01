function Plot_correlation_for_all_chromosomes()
    chr_size = 22;
    K_1_BASE_PATH = '../DATA/Repli_BS/K_RATES/1/';
    K_41_BASE_PATH = '../DATA/Repli_BS/K_RATES/41/';
   OFM = 0.5; %order_of_magnitude
   OUT_FIG_DIR = strcat('Figures/REPLICATE_K_CONSISTENCY_EXP_OFM_',num2str(round(OFM,1)),'/');
   if ~exist(OUT_FIG_DIR)
        mkdir(OUT_FIG_DIR);
   end
   for i = 1 : chr_size
       K1_PATH = strcat(K_1_BASE_PATH, 'chr', num2str(i) , '.mat');
       K41_PATH = strcat(K_41_BASE_PATH, 'chr', num2str(i) , '.mat');
       Figure_path = strcat(OUT_FIG_DIR, 'chr', num2str(i) , '.pdf');
       heatmap(K1_PATH, K41_PATH, Figure_path, OFM, 'log10(K) comb1', 'log10(K) comb41');
   end
end