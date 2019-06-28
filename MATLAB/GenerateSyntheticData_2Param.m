load MLRates_Protocol1a_RandomT_Chr1WT;

NTot=size(FitSites,1); %total number of sites in dataset
NSites=60000; %number of sites to simulate

ChooseSites=datasample(1:NTot,NSites,'Replace',false);
fs=MLEFrac(ChooseSites,1); %assigning f values
ks=MLELam(ChooseSites,1); %assigning k values



ts=[0.5,1.5,4.5,16.5];
for ii=1:numel(ts)
Ps(ii,:)=fs.*(1-exp(-ks*ts(ii))); %probability array 
end

%set the number of reads per site per timepoint. Possibly change this to
%sampling from true read depths
S=NSites;
R0=randi([5,10],S,1);
R1=randi([1,4],S,1);
R2=randi([1,4],S,1);
R3=randi([1,4],S,1);

Numreads=[R0,R1,R2,R3];

%Computing the Read data by probability
for ii=1:S
    for jj=1:numel(ts)
        MReads=sum(double(rand(Numreads(ii,jj),1)<Ps(jj,ii)));
       MethReads(ii,jj)=MReads;
    end
end
UMethReads=Numreads-MethReads;
Times=[0.5,1.5,4.5,16.5];
%next: load MReads and UMethReads into AllDat