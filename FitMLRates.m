function [] = FitMLRates()
DatasetName='First';
datafilename=['AllDat_' DatasetName];
load(datafilename,'AllDat','sites');

tshift=0.5;
Times=[0,1,4,16]+tshift;

LB=0; %lower bound for constrained fit
UB=1000; %upper bound for constrained fit

%AllDat is an array of size (NSites,NTimepoints,2). The data in
%AllDat(i,j,1) is the number of methylated reads at site i at timepoint j.
%The data in AllDat(i,j,2) is the number of unmethylated reads at site i at
%timepoint j. Therefore AllDat(i,j,1)+AllDat(i,j,2)=Totalreads(i,j).

for ii=1:size(AllDat,1); 
    Meths=AllDat(ii,:,1);
    UMeths=AllDat(ii,:,2);
    x0=rand; %it should be checked that results are not sensitive to this initial "guess" parameter
    [x_fmin,f_fmin]=myfmincon(x0);
    
    MLRate(ii)=x_fmin;
end

ratefilename=['MLRates_' DatasetName];
save(ratefilename,'MLRate')

    function [x,fval] = myfmincon(x0)
        fun_fmin=@(x) negLL_fmin(x);
        options=optimset('Display','off');
        [x,fval,exitflag]=fmincon(fun_fmin,x0,[],[],[],[],LB,UB,[],options);
    end

    function dum=negLL_fmin(lambda)
        Pmeth=Meths.*log((1-exp(-lambda.*Times)));
        Pumeth=UMeths.*log(exp(-lambda.*Times));
        LogLikelihood=-sum(Pmeth+Pumeth);
        dum=LogLikelihood;
    end
end




