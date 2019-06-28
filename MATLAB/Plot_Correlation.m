DATA_1_PATH = '../DATA/Repli_BS/K_RATES/1/chr1.mat';%'DATA/k1_chr1.mat';
DATA_41_PATH ='../DATA/Repli_BS/K_RATES/41/chr1.mat'; %'DATA/k41_chr1.mat';
    
load(DATA_1_PATH);
k1 = zeros(length(FitSites'), 2); 
k1(:, 1) = FitSites';% Fitted sites
k1(:, 2) = MLELam(:, 1); % K values

clear FitSites MLELam MLEFrac;

load(DATA_41_PATH);

k41 = zeros(length(FitSites'), 2);
k41(:, 1) = FitSites';
k41(:, 2) = MLELam(:, 1);

clear FitSites MLELam MLEFrac;

[overlapped_sites, ia, ib] = intersect(k1(:, 1), k41(:, 1)); % # of Sites in experiment 1 and 41 are: 153191 and 160712, #overlapped site is 46038 ~30% of exp1, 28% of exp41.

overlapped_k1 = log(k1(ia, 2)); % Log of the K values of the overlapped sites in data 1.
overlapped_k41 = log(k41(ib, 2)); % in data 41

over_idxs = overlapped_k1 >= -2 & overlapped_k41 >= -2; % Get rid of the -Inf values in the logged K for both data 1 and data 41.
overlapped_k1=overlapped_k1(over_idxs); % -Inf filtered overlapped K values in data 1
overlapped_k41=overlapped_k41(over_idxs);% -Inf filtered overlapped K values in data 41

fig = figure(1);

subplot(2,2,1);
X_edges = -2: 0.1: 1; % The range of K values in data 1
Y_edges = -2: 0.1: 1; % The range of K values in data 41

[N,C]=hist3([overlapped_k1(:), overlapped_k41(:)], 'Edges', {X_edges, Y_edges});
contourf(C{1},C{2},N); % Contour Plot of K values in both data 1 and data 41
colorbar;

[R,P,RLO,RUP]= corrcoef(overlapped_k1, overlapped_k41, 'alpha', 0.05); % Correlation coefficient of K values between data 1 and data 41


tt1 = strcat('Corrlation: ', num2str(round(R(2,1),2)), ', 95% CI: [', num2str(round(RLO(2,1),2)), ',' , num2str(round(RUP(2,1),2)) , ']');
title(tt1);
xlabel('K from DATA 1');
ylabel('K from DATA 41');



subplot(2,2,2);
close_idxs = abs(overlapped_k1 - overlapped_k41) < 1; % Whether the K values in overlapped CpGs are close enough (within 1 order of mangnitude, i.e. abs(Ks in 1 - Ks in 41) < 1), return a boolean array
far_idxs = ~close_idxs; % The boolean array for the CpGs that are not close enough, return the indexs of them.

percentage_of_close_idxs = double(length(overlapped_k1(close_idxs)))/length(overlapped_k41) * 100;

plot(overlapped_k1(close_idxs), overlapped_k41(close_idxs), 'r.');
[R2,P2,RLO2,RUP2]= corrcoef(overlapped_k1(close_idxs), overlapped_k41(close_idxs), 'alpha', 0.05);
hold on;
plot(overlapped_k1(far_idxs), overlapped_k41(far_idxs), 'b.');
tt = strcat(num2str(length(overlapped_k1(close_idxs))) ,'# , ', num2str(round(percentage_of_close_idxs, 1)), '%  in red, Corr: ', num2str(round(R2(2,1),2)), ', 95% CI: [', num2str(round(RLO2(2,1),2)), ',' , num2str(round(RUP2(2,1),2)) , ']');

title(tt);
xlabel('K from DATA 1');
ylabel('K from DATA 41');
xlim([-2, 1]);
ylim([-2, 1]);

FIGURE_DIR = 'Figures';
if ~exist(FIGURE_DIR)
    mkdir(FIGURE_DIR);
end
print(fig, strcat(FIGURE_DIR, '/', 'Density_of_logK_in_data_1_and_41.pdf') , '-dpdf','-opengl','-r300');
close;