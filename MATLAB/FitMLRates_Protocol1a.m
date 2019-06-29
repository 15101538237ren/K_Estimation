function [] = FitMLRates_Protocol1a(data_fp, out_fp)
%read in the dataset
load(data_fp, 'AllDat', 'sites');

%the experimental timepoints. They are shifted by tshift, which accounts
%for the fact that at, e.g., time=0 hrs, the captured reads started
%replicating within the previous hour. This is done also because the model
%assumes zero probability of methylation occuring by time=0.
tshift=0.5;
Times=[0,1,4,16]+tshift;

%Select only the sites which have 'enough' data--a certain number of
%collected reads and timepoints
EnoughReads0=10; %number of reads required at t=0
EnoughReadsLater=5; %number of reads required for later timepoints
Reads=sum(AllDat(:,:,1:2),3);
NumReads0=Reads(:,1);
NumReadsLater=sum(Reads(:,2:end),2);
KeepSites=find(NumReadsLater>=EnoughReadsLater & NumReads0>=EnoughReads0);
numsites=numel(KeepSites);

logminlam=-2; %log10 minimum allowed fit rate
logmaxlam=1;  %log10 maximume allowed fit rate--calculated according to ReadDepth
fracgrid=[0:0.01:1];    %grid of f-values at which LogLikelihood will be computed
lamgrid=10.^([logminlam:.01:logmaxlam]);   %grid of k-values
maxlam=max(lamgrid);
minlam=min(lamgrid);

%initialize various arrays
MLELam=zeros(numsites,3);
MLEFrac=zeros(numsites,3);
FitSites=sites(KeepSites);
Fitness=zeros(numsites,1);

tic
SaveEvery=5E4;
counter=0;

parfor ii=1:numsites %loop over the sites
    
    PmethA=zeros(numel(fracgrid),numel(lamgrid),numel(Times));
    PumethA=zeros(numel(fracgrid),numel(lamgrid),numel(Times));
    
    %save the fitted data periodically
%     if (ii-counter*SaveEvery)/SaveEvery >= 1
%         save ii ii
%         counter=counter+1;
%         save(out_fp,'FitSites','MLELam','MLEFrac')
%         toc
%     end

    
    %get the read data for this site
    Meths=AllDat(KeepSites(ii),:,1);
    UMeths=AllDat(KeepSites(ii),:,2);
    
    %loop over the timepoints to compute the LogLikelihood surface
    for tind=1:numel(Times)
        Time=Times(tind);
%for time tj and parameters f and k (lambda), the probability of observing methyl is assumed to be
%P=f+f/k*(exp(-k*(tj+.5))-exp(-k*(tj-0.5)))
        expterms=1-exp(-lamgrid.*Time);
        expterms(expterms>=1)=1-eps;
        repexpterms=repmat(expterms,numel(fracgrid),1);
        repfrac=repmat(fracgrid',1,numel(lamgrid));
        PmethA(:,:,tind)=Meths(tind).*log(repfrac.*repexpterms);
        PumethA(:,:,tind)=UMeths(tind).*log(1-repfrac.*repexpterms);
    end
    %this is the LogLikelihood surface as a function of the parameters
    LogLikelihood=sum(PmethA,3,'omitnan')+sum(PumethA,3,'omitnan');
    
    %find the parameter values that maximize the LogLikelihood
    szL=size(LogLikelihood);
    [MaxLL,i]=max(LogLikelihood(:));
    [I,J]=ind2sub(szL,i); %I,J are the indices of ML parameters
    %these are the "raw" ML parameters, but need to deal with edge cases (below)	
    RawLam=lamgrid(J);
    RawFrac=fracgrid(I);
    
    %compute the profile likelihood functions for both parameters. (These are used
    %to compute confidence intervals).
    ProfileFrac=max(LogLikelihood,[],2);
    ProfileLam=max(LogLikelihood,[],1);
    GetCILam=-2*(ProfileLam-LogLikelihood(I,J));
    GetCIfrac=-2*(ProfileFrac-LogLikelihood(I,J));
    
    %CIs will be computed based on Likelihood Ratio Test, using Chi-squared
    %values
    PrcVal=3.841; %95th of chi-squared distribution, 1 free param. 
    PrcValIn=1.32; %75th
    
    %Find the 95% confidence intervals
    LamRegion95=lamgrid(GetCILam<=PrcVal);
    fracRegion95=fracgrid(GetCIfrac<=PrcVal);
    CILam95=[LamRegion95(1) LamRegion95(end)];
    CIfrac95=[fracRegion95(1) fracRegion95(end)];
    
    %Find the inner confidence intervals (as determined by PrcValIn)
    LamRegionIn=lamgrid(GetCILam<=PrcValIn);
    fracRegionIn=fracgrid(GetCIfrac<=PrcValIn);
    CILamIn=[LamRegionIn(1) LamRegionIn(end)];
    CIfracIn=[fracRegionIn(1) fracRegionIn(end)];

    %Find the indices in lamgrid of the inner CI edges (upper and lower)
    LamIndsIn=find(GetCILam<=PrcValIn);
    LamIndsInnerEdges=[LamIndsIn(1) LamIndsIn(end)];
    
    %deal with edge cases: where k is either unidentifiable, or only lower/upper bound is identifiable
    %check to see whether the inner confidence interval hits the boundary
    CheckCIsL=[minlam,maxlam]-CILamIn;
    keepindlam=find(abs(CheckCIsL)>realmin);
    
    if numel(keepindlam)==0 %this indicates the inner CI overlaps entire region-->k is unidentifiable. Set k==0
        %this was tested--the sites with only unmethylated reads or very few methylated reads are found
        %this way
        MLELam(ii,:)=[0,0,0];
        MLEFrac(ii,:)=[0,0,0];
    elseif abs(RawLam-maxlam)<=realmin %if the raw k lies on upper boundary of range
        %set k value to the lower edge of inner CI (gives a lower bound on k)
	MLELam(ii,:)=[CILamIn(1),CILam95(1),CILam95(2)];
        lamind=LamIndsInnerEdges(1);
        %now find the corresponding max f value
        [aa,bb]=max(LogLikelihood(:,lamind));
        MLEFrac(ii,:)=[fracgrid(bb),CIfrac95(1),CIfrac95(2)];
        MaxLL=LogLikelihood(bb,lamind);
    elseif abs(RawLam-minlam)<=realmin %if the raw k lies on lower boundary of range
        %set k value to the upper edge of inner CI (gives an upper bound on k)
	MLELam(ii,:)=[CILamIn(2),CILam95(1),CILam95(2)];
        lamind=LamIndsInnerEdges(2);
        %now find the corresponding max f value
        [aa,bb]=max(LogLikelihood(:,lamind));
        MLEFrac(ii,:)=[fracgrid(bb),CIfrac95(1),CIfrac95(2)];
        MaxLL=LogLikelihood(bb,lamind);
    else %if no edge case is relevant
        MLELam(ii,:)=[RawLam,CILam95(1),CILam95(2)];
        MLEFrac(ii,:)=[RawFrac,CIfrac95(1),CIfrac95(2)];
    end
Fitness(ii)=MaxLL;
end
save(out_fp,'FitSites','MLELam','MLEFrac')
toc
end



