clear;
clc;
close all;
K_1_PATH = '../DATA/Repli_BS/K_RATES/1/chr1.mat';%'DATA/k1_chr1.mat';
load(K_1_PATH);
NTot=size(FitSites, 2); %total number of sites in dataset
NSites=60000; %number of sites to simulate

rs = RandStream('mlfg6331_64');%Create the random number stream for reproducibility.

ChooseSites=datasample(rs, 1 : NTot, NSites, 'Replace',false);
fs=MLEFrac(ChooseSites,1); %assigning f values
ks=MLELam(ChooseSites ,1); %assigning k values

ts=[0.5,1.5,4.5,16.5];
N_Times= numel(ts);
for ii=1: N_Times
Ps(ii,:)=fs.*(1-exp(-ks*ts(ii))); %probability array 
end


MERGED_CHR1_DATA_PATH = '../DATA/Repli_BS/TMP/MERGED_DATA/chr1.mat';
load(MERGED_CHR1_DATA_PATH, 'AllDat', 'sites'); % All merged chr1 sites, without Fitting, therefore, to sample the reads from the fitted sites, we have to extract the fitted sites and sample from them.
[olps, ia, ib] = intersect(FitSites, sites); % !!!!!!!!This probably not appropriate.
Fitted_dat_with_reads_info = sum(AllDat(ib, :, 1:2), 3);

%set the number of reads per site per timepoint. Possibly change this to
%sampling from true read depths
Tot_Size = size(Fitted_dat_with_reads_info, 1);
sampled_idxs = datasample(rs, 1 : Tot_Size, NSites, 'Replace',false);
Sampled_data_with_reads_info = Fitted_dat_with_reads_info(sampled_idxs , :);


All_Reads_rep1 = zeros(NSites, N_Times, 2);
All_Reads_rep2 = zeros(NSites, N_Times, 2);

p = 0.5; %probability of binomial distribution
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

REP1_DATA_PATH = 'DATA/REP1.mat';
AllDat = All_Reads_rep1;
save(REP1_DATA_PATH,'AllDat','sites');

REP1_K_PATH = 'DATA/REP1_K.mat';
FitMLRates_Protocol1a_RandomT(REP1_DATA_PATH, REP1_K_PATH);

REP2_DATA_PATH = 'DATA/REP2.mat';
AllDat = All_Reads_rep2;
save(REP2_DATA_PATH,'AllDat','sites');

REP2_K_PATH = 'DATA/REP2_K.mat';
FitMLRates_Protocol1a_RandomT(REP2_DATA_PATH, REP2_K_PATH);