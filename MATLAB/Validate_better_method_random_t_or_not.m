clear;
clc;
close all;

OFM = 0.5; % The order of magnitude for which K difference can be considered as conservative.

chr_size = 1;
NSites=60000; %number of sites to simulate
p = 0.5; %probability of binomial distribution
ts=[0.5,1.5,4.5,16.5];
N_Times= numel(ts);
random_ts = rand(NSites, 4) + ts - 0.5; % Add random T when generates the data
rs = RandStream('mlfg6331_64');%Create the random number stream for reproducibility.

K_1_BASE_PATH = '../DATA/Repli_BS/K_RATES/1/';
MERGED_DATA_PATH = '../DATA/Repli_BS/TMP/MERGED_DATA/';
OUT_FIG_DIR = strcat('Figures/INFER_METHOD_COMPARISON_OFM_',num2str(round(OFM,1)),'/');
if ~exist(OUT_FIG_DIR)
    mkdir(OUT_FIG_DIR);
end
for i = 1 : chr_size
    K_1_PATH = strcat(K_1_BASE_PATH, 'chr', num2str(i) , '.mat');
    load(K_1_PATH);
    NTot=size(FitSites, 2); %total number of sites in dataset
    
    ChooseSites=datasample(rs, 1 : NTot, NSites, 'Replace',false);
    F_GROUND_TRUTH= MLEFrac(ChooseSites,:);
    K_GROUND_TRUTH= MLELam(ChooseSites,:);
    fs=MLEFrac(ChooseSites,1); %assigning f values
    ks=MLELam(ChooseSites ,1); %assigning k values

    Ps = zeros(N_Times, NSites);
    for ii=1: N_Times
        Ps(ii,:)=fs.*(1-exp(-ks.*random_ts(:, ii))); %probability array 
    end
    
    load(strcat(MERGED_DATA_PATH, 'chr', num2str(i) , '.mat'), 'AllDat', 'sites'); % All merged chr1 sites, without Fitting, therefore, to sample the reads from the fitted sites, we have to extract the fitted sites and sample from them.

    % The following two lines selected the reads only by the fitted sites
    [olps, ia, ib] = intersect(FitSites, sites);
    Fitted_dat_with_reads_info = sum(AllDat(ib, :, 1:2), 3);
    origin_data_with_reads_info = sum(AllDat(:, :, 1:2), 3);
    MEAN_DEPTH_AT_DIFF_POINTS_origin = mean(origin_data_with_reads_info / 2.0, 1);
    STD_DEPTH_AT_DIFF_POINTS_origin = std(origin_data_with_reads_info / 2.0, 1);
    

    %Fitted_dat_with_reads_info = sum(AllDat(:, :, 1:2), 3);

    %set the number of reads per site per timepoint. Possibly change this to
    %sampling from true read depths
    Tot_Size = size(Fitted_dat_with_reads_info, 1);
    sampled_idxs = datasample(rs, 1 : Tot_Size, NSites, 'Replace',false);
    Sampled_data_with_reads_info = Fitted_dat_with_reads_info(sampled_idxs , :);


    All_Reads_rep1 = zeros(NSites, N_Times, 2);
    All_Reads_rep2 = zeros(NSites, N_Times, 2);

    
    %Computing the Read data by probability
    for ii=1 : NSites
        for jj=1 : N_Times
            Numreads = Sampled_data_with_reads_info(ii,jj);
            MReads=sum(double(rand(Numreads, 1) < Ps(jj,ii)));
            UReads=Numreads-MReads;
            MR1 = binornd(MReads, p);
            All_Reads_rep1(ii,jj, 1)=MR1;
            All_Reads_rep2(ii,jj, 1)=MReads - MR1;

            UR1 = binornd(UReads, p);
            All_Reads_rep1(ii,jj, 2)=UR1; % Unmethylated Reads rep1
            All_Reads_rep2(ii,jj, 2)=UReads - UR1;
        end
    end
    sites = FitSites(ChooseSites);

    REP1_DATA_PATH = 'DATA/REP1.mat'; % Save mat data of rep1 and rep2
    AllDat = All_Reads_rep1;
    save(REP1_DATA_PATH,'AllDat','sites');
    
    REP2_DATA_PATH = 'DATA/REP2.mat';
    AllDat = All_Reads_rep2;
    save(REP2_DATA_PATH,'AllDat','sites');
    
    REP1_K_NON_RANDOM = 'DATA/REP1_K_NON_RANDOM.mat';
    REP1_K_RANDOM = 'DATA/REP1_K_RANDOMT.mat';
    
    REP2_K_NON_RANDOM = 'DATA/REP2_K_NON_RANDOM.mat';
    REP2_K_RANDOM = 'DATA/REP2_K_RANDOMT.mat';
    
    %Infer K of rep1 by non-random t random t also for rep2
%      FitMLRates_Protocol1a(REP1_DATA_PATH, REP1_K_NON_RANDOM);
%      FitMLRates_Protocol1a_RandomT(REP1_DATA_PATH, REP1_K_RANDOM);
%      FitMLRates_Protocol1a(REP2_DATA_PATH, REP2_K_NON_RANDOM);
%      FitMLRates_Protocol1a_RandomT(REP2_DATA_PATH, REP2_K_RANDOM);
    
    MLELam = K_GROUND_TRUTH;
    MLEFrac = F_GROUND_TRUTH;
    FitSites = sites;
    GROUND_TRUTH_K_fp = 'DATA/GROUND_TRUTH_K.mat';
    save(GROUND_TRUTH_K_fp,'FitSites','MLELam','MLEFrac');
    
    label1 = 'log10(k) ground truth';
    random_label1 = 'log10(k) rep1 random t';
    random_label2 = 'log10(k) rep2 random t';
    non_random_label1 = 'log10(k) rep1 non-random t';
    non_random_label2 = 'log10(k) rep2 non-random t';
%      REP1_K_NON_RANDOM_fig = strcat(OUT_FIG_DIR, 'REP1_K_NON_RANDOM.pdf');
%      heatmap(GROUND_TRUTH_K_fp, REP1_K_NON_RANDOM, REP1_K_NON_RANDOM_fig, OFM, label1, non_random_label1);
% %     
%      REP2_K_NON_RANDOM_fig = strcat(OUT_FIG_DIR, 'REP2_K_NON_RANDOM.pdf');
%      heatmap(GROUND_TRUTH_K_fp, REP2_K_NON_RANDOM, REP2_K_NON_RANDOM_fig, OFM, label1, non_random_label2);
% %     
     REP1_K_RANDOM_fig = strcat(OUT_FIG_DIR, 'REP1_K_RANDOM.pdf');
     heatmap(GROUND_TRUTH_K_fp, REP1_K_RANDOM, REP1_K_RANDOM_fig, OFM, label1, random_label1);
%     
     REP2_K_RANDOM_fig = strcat(OUT_FIG_DIR, 'REP2_K_RANDOM.pdf');
     heatmap(GROUND_TRUTH_K_fp, REP2_K_RANDOM, REP2_K_RANDOM_fig, OFM, label1, random_label2);
%     
%     
%     REP1_K_RANDOM_DEPTH_fig = strcat(OUT_FIG_DIR, 'REP1_K_RANDOM_DEPTH.pdf');
%     plot_k_estimation_with_read_depth(GROUND_TRUTH_K_fp, REP1_K_RANDOM, REP1_DATA_PATH, REP1_K_RANDOM_DEPTH_fig, label1, random_label1);
%     
%     REP2_K_RANDOM_DEPTH_fig = strcat(OUT_FIG_DIR, 'REP2_K_RANDOM_DEPTH.pdf');
%     plot_k_estimation_with_read_depth(GROUND_TRUTH_K_fp, REP2_K_RANDOM, REP2_DATA_PATH, REP2_K_RANDOM_DEPTH_fig, label1, random_label2);
%     
    %Improve the read depth in synthesized data by adding constant read
    %depth D
%     MEAN_DEPTH_AT_DIFF_POINTS = mean(Sampled_data_with_reads_info / 2.0, 1);
%     STD_DEPTH_AT_DIFF_POINTS = std(Sampled_data_with_reads_info / 2.0, 1);
%     READ_DEPTH_FIG_DIR = strcat(OUT_FIG_DIR, 'DIFF_READ_DEPTH_INCRESE_AT_1h_4h_16h/');
%     if ~exist(READ_DEPTH_FIG_DIR)
%         mkdir(READ_DEPTH_FIG_DIR);
%     end
%     for D = 5: 5: 20
%         All_Reads_rep1 = zeros(NSites, N_Times, 2);
%         All_Reads_rep2 = zeros(NSites, N_Times, 2);
%         for ii=1 : NSites
%             for jj=1 : N_Times
%                 Numreads = Sampled_data_with_reads_info(ii,jj);
%                 if jj ~= 1
%                      Numreads = Numreads + D * 4;
%                 end
%                 MReads=sum(double(rand(Numreads, 1) < Ps(jj,ii)));
%                 UReads=Numreads-MReads;
%                 MR1 = binornd(MReads, p);
%                 All_Reads_rep1(ii,jj, 1)=MR1;
%                 All_Reads_rep2(ii,jj, 1)=MReads - MR1;
% 
%                 UR1 = binornd(UReads, p);
%                 All_Reads_rep1(ii,jj, 2)=UR1; % Unmethylated Reads rep1
%                 All_Reads_rep2(ii,jj, 2)=UReads - UR1;
%             end
%         end
%         REP1_DATA_PATH = 'DATA/REP1.mat'; % Save mat data of rep1 and rep2
%         AllDat = All_Reads_rep1;
%         save(REP1_DATA_PATH,'AllDat','sites');
% 
%         REP2_DATA_PATH = 'DATA/REP2.mat';
%         AllDat = All_Reads_rep2;
%         save(REP2_DATA_PATH,'AllDat','sites');
%         
%         REP1_K_RANDOM = strcat('DATA/REP1_K_RANDOMT_',num2str(D),'.mat');
%          FitMLRates_Protocol1a_RandomT(REP1_DATA_PATH, REP1_K_RANDOM);
%         REP2_K_RANDOM = strcat('DATA/REP2_K_RANDOMT_',num2str(D),'.mat');
%          FitMLRates_Protocol1a_RandomT(REP2_DATA_PATH, REP2_K_RANDOM);
%         
%         REP1_K_RANDOM_fig = strcat(READ_DEPTH_FIG_DIR, 'REP1_K_RANDOM_',num2str(D),'.pdf');
%         heatmap(GROUND_TRUTH_K_fp, REP1_K_RANDOM, REP1_K_RANDOM_fig, OFM, label1, random_label1);
%   
%         REP2_K_RANDOM_fig = strcat(READ_DEPTH_FIG_DIR, 'REP2_K_RANDOM_',num2str(D),'.pdf');
%         heatmap(GROUND_TRUTH_K_fp, REP2_K_RANDOM, REP2_K_RANDOM_fig, OFM, label1, random_label2);
%     end
end

