function [] = FitMLRates_Protocol1()
%read in the dataset
DatasetName=['Sim_DiffDec_WTBulkf_D2'];
ProtocolName='Protocol1';
datafilename=['AllDat_' DatasetName];
load(datafilename,'AllDat','sites');
%AllDat is an array of size (NSites,NTimepoints,2). The data in
%AllDat(i,j,1) is the number of methylated reads at site i at timepoint j.
%The data in AllDat(i,j,2) is the number of unmethylated reads at site i at
%timepoint j. Therefore AllDat(i,j,1)+AllDat(i,j,2)=Totalreads(i,j).

ratefilename=['MLRates_' ProtocolName '_' DatasetName]; %filename for saving results

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
FitSites=zeros(numsites,1);
PmethA=zeros(numel(fracgrid),numel(lamgrid),numel(Times));
PumethA=zeros(numel(fracgrid),numel(lamgrid),numel(Times));

tic
SaveEvery=5E4;
counter=0;

for ii=1:numsites %loop over the sites
    
    %save the fitted data periodically
    if (ii-counter*SaveEvery)/SaveEvery >= 1
        save ii ii
        counter=counter+1;
        save(ratefilename,'FitSites','MLELam','MLEFrac')
        toc
    end
    
    %get the read data for this site
    Meths=AllDat(KeepSites(ii),:,1);
    UMeths=AllDat(KeepSites(ii),:,2);
    
    %loop over the timepoints to compute the LogLikelihood surface
    for tind=1:numel(Times)
        Time=Times(tind);
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
    [I,J]=ind2sub(szL,i);
    
    %compute the profile likelihood functions for both parameters. (These are used
    %to compute confidence intervals).
    ProfileFrac=max(LogLikelihood,[],2);
    ProfileLam=max(LogLikelihood,[],1);
    GetCILam=-2*(ProfileLam-LogLikelihood(I,J));
    GetCIfrac=-2*(ProfileFrac-LogLikelihood(I,J));
    
    %CIs will be computed based on Likelihood Ratio Test, using Chi-squared
    %values
    PrcVal=3.841; %95th of chi-squared distribution, 1 free param. 
    %PrcValIn=PrcVal;%95th of chi-squared distribution, 1 free param.
    PrcValIn=1.32; %75th
    %PrcValIn=0.455; %50th
    %PrcValIn=0.102; %25th
    %PrcValIn=0.016; %10th
    
    %Find the 95% confidence intervals
    LamRegion95=lamgrid(GetCILam<=PrcVal);
    fracRegion95=fracgrid(GetCIfrac<=PrcVal);
    CILam95=[LamRegion95(1) LamRegion95(end)];
    CIfrac95=[fracRegion95(1) fracRegion95(end)];
    
    %Find the inner (50th%) confidence intervals
    LamRegionIn=lamgrid(GetCILam<=PrcValIn);
    fracRegionIn=fracgrid(GetCIfrac<=PrcValIn);
    CILamIn=[LamRegionIn(1) LamRegionIn(end)];
    CIfracIn=[fracRegionIn(1) fracRegionIn(end)];
    
    %check to see whether the inner confidence interval hits the boundary
    CheckCIsL=[minlam,maxlam]-CILamIn;
    keepindlam=find(abs(CheckCIsL)>realmin);
    
    if numel(keepindlam)==0 %this indicates the inner CI overlaps entire region-->k is unidentifiable. Set k==0
        %this was tested--the sites with only unmethylated reads are found
        %this way, and this method gives identical results to simply
        %searching for reads with only unmethylated sites
        MLELam(ii,:)=[0,0,0];
        MLEFrac(ii,:)=[0,0,0];
    elseif numel(keepindlam)==1 %this indicates CI value hits one boundary-->indicated upper/lower bound on k
        %when the CI hits the boundary, use the inner CI value as the rate
        %estimate
        MLELam(ii,:)=[CILamIn(keepindlam),CILam95(1),CILam95(2)];
        lamind=find(lamgrid==CILamIn(keepindlam));
        [aa,bb]=max(LogLikelihood(:,lamind));
        MLEFrac(ii,:)=[fracgrid(bb),CIfrac95(1),CIfrac95(2)];
    else
        MLELam(ii,:)=[lamgrid(J),CILam95(1),CILam95(2)];
        MLEFrac(ii,:)=[fracgrid(I),CIfrac95(1),CIfrac95(2)];
    end
    FitSites(ii)=sites(KeepSites(ii));
end
save(ratefilename,'FitSites','MLELam','MLEFrac')
toc
end




