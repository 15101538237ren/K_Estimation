function Threshold_impact_on_number_of_overlapped_CpGs()
    DATA_1_PATH = '../DATA/Repli_BS/MATLAB_DATA/1/chr1.mat';%'DATA/origin_1_chr1.mat';
    DATA_41_PATH = '../DATA/Repli_BS/MATLAB_DATA/41/chr1.mat'; %'DATA/origin_41_chr1.mat';
    thresholds_1 = 5 : 15; % Thresholds for the number of reads required at t=0
    thresholds_2 = 1 : 7; % Thresholds for the number of reads required at later timepoints
    
    l1 = length(thresholds_1);
    l2 = length(thresholds_2);
    
    load(DATA_1_PATH, 'AllDat', 'sites');
    Reads_of_1=sum(AllDat(:,:,1:2),3);
    Sites_of_1 = sites;
    clear AllDat sites;
    
    load(DATA_41_PATH, 'AllDat', 'sites');
    Reads_of_41=sum(AllDat(:,:,1:2),3);
    Sites_of_41 = sites;
    clear AllDat sites;
    
    number_of_total_sites = double(length(Sites_of_1));
    
    percentages_1 = zeros(l1, l2); % The percentage matrix of the CpGs Filtered by threshold but not overlapped in data 1.
    percentages_41 = zeros(l1, l2); % The percentage matrix of the CpGs Filtered by threshold but not overlapped in data 41.
    percentages_overlapped = zeros(l1, l2); % Initialize the percentage matrix of the overlapped CpGs with different thresholds
    ovlp_count_10_5 = 0; % The overlapped CpGs for threshold 1 = 10 and threshold2 = 5.
    for i = 1: l1
        for j = 1 : l2
            th1 = thresholds_1(i);
            th2 = thresholds_2(j);
            
            sites_1 = Filter_the_sites_by_thresholds(Reads_of_1, Sites_of_1, th1, th2);
            sites_41 = Filter_the_sites_by_thresholds(Reads_of_41, Sites_of_41, th1, th2);
            
            percentages_1(i, j) = length(sites_1) / number_of_total_sites;
            percentages_41(i, j) = length(sites_41) / number_of_total_sites;
            
            overlapped_sites = intersect(sites_1, sites_41); % The overlapped sites between data 1 and data 41
            if th1 == 10 && th2 ==5
                ovlp_count_10_5 = length(overlapped_sites);
            end
            percentages_overlapped(i, j) = length(overlapped_sites) / number_of_total_sites;
        end
    end
    ratio_of_overlapped_cpg = ovlp_count_10_5/number_of_total_sites;
    fig = figure(1);
    subplot(2,2,1);
    imagesc(percentages_overlapped);
    configuration_for_each_subplot();
    title(strcat('%Overlapped CpGs: (10, 5)= ', num2str(ovlp_count_10_5), '/', num2str(number_of_total_sites), '=', num2str(round(ratio_of_overlapped_cpg,2))));
    
    subplot(2,2,2);
    imagesc(percentages_1);
    configuration_for_each_subplot();
    title('%CpGs in 1 Filtered, No Overlap');
    
    subplot(2,2,3);
    imagesc(percentages_41);
    configuration_for_each_subplot();
    title('%CpGs in 41 Filtered, No Overlap');
    
    FIGURE_DIR = 'Figures';
    if ~exist(FIGURE_DIR)
        mkdir(FIGURE_DIR);
    end
    print(fig, strcat(FIGURE_DIR, '/', 'Threshold.pdf'), '-dpdf','-opengl','-r300');
    close;
end