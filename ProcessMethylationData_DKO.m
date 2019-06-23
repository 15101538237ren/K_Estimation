%function [AllDat,PercMeth] = ProcessMethylationData_DKO()
filename='DKO/DKO.xlsx';
DatasetName='DKO';
datafilename=['AllDat_' DatasetName];%'TiledDat_500';
data=xlsread(filename,1);
sites=data(:,1);
for ii=2:7
data=xlsread(filename,ii);
sites=unique([sites;data(:,1)]); %just list all sites
end
NSites=numel(sites);
NTimePoints=4;
AllDat=zeros(NSites,NTimePoints,2);
Times=[1,1,2,2,3,3,4,4];
for ii=1:7
	data=xlsread(filename,ii);
NMeth=data(:,4);
NUMeth=data(:,5)-NMeth;
[Lia,Locb]=ismember(data(:,1),sites);
AllDat(Locb,Times(ii),1)=AllDat(Locb,Times(ii),1)+NMeth;
AllDat(Locb,Times(ii),2)=AllDat(Locb,Times(ii),2)+NUMeth;
end

%AllDat_DKO=AllDat;
%sites_DKO=sites;

save(datafilename,'AllDat','sites');
%sites_DKO sites_DKO
%save AllDat_DKO AllDat_DKO



