function heatmap(fitfile1, fitfile2, figure_path, order_of_magnitude)
    %this script looks at correlation between fitted rates (k) from different
    %methods

    load(fitfile1,'FitSites','MLELam');

    Sites1=FitSites;
    Rates1=MLELam(:,1);
    load(fitfile2,'FitSites','MLELam');

    Sites2=FitSites;
    Rates2=MLELam(:,1);%FitLam;
    
%     fig2 = figure(2);
%     subplot(2,2,1);
%     histogram(log10(Rates1));
%     subplot(2,2,2);
%     histogram(log10(Rates2));

    % %find the sites that intersect (should be all of them)
    [C,ia,ib]=intersect(Sites1,Sites2);
    %%Build up a 2D histogram comparing the two datasets
    logminlam=-2; %log10 minimum allowed fit rate
    logmaxlam=1;  %log10 maximume allowed fit rate--calculated according to ReadDepth
    lamgrid=10.^([logminlam:.1:logmaxlam]);   %grid of k-values
    LamArray=zeros(numel(lamgrid));
    for ii=1:numel(C)
        indsa=find(Rates1(ia(ii))<=lamgrid);
        indsb=find(Rates2(ib(ii))<=lamgrid);
        LamArray(indsa(1),indsb(1))=LamArray(indsa(1),indsb(1))+1;
    end
    
    LamArray=LamArray/numel(C);
    CC=corrcoef(Rates1(ia),Rates2(ib));
    lala=log10(lamgrid);
    fig = figure(1);
    
    subplot(2,2,1);
    imagesc(lala,lala,-log(LamArray));
    colormap(flipud(parula));
    h=colorbar;
    ylabel(h,'Probability')
    Vals=fliplr([0.0001,0.001,0.01,0.1]);
    Ticks=-log(Vals);
    set(h,'YDir','reverse');
    set(h,'Ticks',Ticks);
    set(h,'TickLabels',cellstr(num2str(Vals(:)))');
    axis square;
    xlabel('log10(k) DATA 1');
    ylabel('log10(k) DATA 2');
    set(gca, 'YDir', 'normal');
    title(strcat('Corrlation: ', num2str(round(CC(2,1),2))));
    
    subplot(2,2,2);
    log10_ovlp_K1 = log10(Rates1(ia)); % Log of the K values of the overlapped sites in data 1.
    log10_ovlp_K2 = log10(Rates2(ib)); % in data 41

    ovlp_idxs = log10_ovlp_K1 >= -2 & log10_ovlp_K2 >= -2; % Get rid of the -Inf values in the logged K for both data 1 and data 41.
    log10_ovlp_K1=log10_ovlp_K1(ovlp_idxs); % -Inf filtered overlapped K values in data 1
    log10_ovlp_K2=log10_ovlp_K2(ovlp_idxs);% -Inf filtered overlapped K values in data 41

    conservative_indexs = abs(log10_ovlp_K1 - log10_ovlp_K2) < order_of_magnitude; % Whether the K values in overlapped CpGs are close enough (within 1 order of mangnitude, i.e. abs(Ks in 1 - Ks in 41) < 1), return a boolean array

    ratio_of_conservative = double(length(log10_ovlp_K1(conservative_indexs)))/length(log10_ovlp_K1) * 100.0;
    
    conservative_K1 = log10_ovlp_K1(conservative_indexs);
    conservative_K2 = log10_ovlp_K2(conservative_indexs);
    unconservative_K1 = log10_ovlp_K1(~conservative_indexs);
    unconservative_K2 = log10_ovlp_K2(~conservative_indexs);
    
    plot(conservative_K1, conservative_K2, 'r.');
    hold on;
    plot(unconservative_K1, unconservative_K2, 'b.');
    
    tt = strcat(num2str(length(conservative_K1)) ,' and ', num2str(round(ratio_of_conservative, 1)), ' % conservative K in overlapped sites');
    title(tt);
    xlabel('log10(k) DATA 1');
    ylabel('log10(k) DATA 2');
    xlim([-2, 1]);
    ylim([-2, 1]);
    print(fig, figure_path, '-dpdf','-opengl','-r300');
    close;
end


