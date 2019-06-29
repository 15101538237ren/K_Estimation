function Plot_correlation_for_all_chromosomes()
    chr_size = 22;
    K_1_BASE_PATH = '../DATA/Repli_BS/K_RATES_NON_RANDOM/1/';
    K_41_BASE_PATH = '../DATA/Repli_BS/K_RATES_NON_RANDOM/41/';
   OUT_FIG_DIR = 'Figures/EXPERIMENT_1_and_41_0.5_ORDER_OF_MAGNITITUDE_NON_RANDOM_K_MLE/';
   order_of_magnitude = 0.5;
   if ~exist(OUT_FIG_DIR)
        mkdir(OUT_FIG_DIR);
   end
   for i = 1 : chr_size
       K1_PATH = strcat(K_1_BASE_PATH, 'chr', num2str(i) , '.mat');
       K41_PATH = strcat(K_41_BASE_PATH, 'chr', num2str(i) , '.mat');
       Figure_path = strcat(OUT_FIG_DIR, 'chr', num2str(i) , '.pdf');
       heatmap(K1_PATH, K41_PATH, Figure_path, order_of_magnitude);
   end
end