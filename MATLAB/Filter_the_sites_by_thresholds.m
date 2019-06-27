function [FitSites] = Filter_the_sites_by_thresholds(Reads, sites, EnoughReads0, EnoughReadsLater)
NumReads0=Reads(:,1);
NumReadsLater=sum(Reads(:,2:end),2);
KeepSites=NumReadsLater>=EnoughReadsLater & NumReads0>=EnoughReads0;
FitSites=sites(KeepSites); % The filtered sites to return
end