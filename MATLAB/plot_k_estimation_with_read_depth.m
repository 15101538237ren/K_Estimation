function plot_k_estimation_with_read_depth(fitfile1, fitfile2, read_file2, figure_path, label1, label2)
    load(fitfile1,'FitSites','MLELam');
    Sites1=FitSites;
    Rates1=MLELam(:,1);

    load(fitfile2,'FitSites','MLELam');
    Sites2=FitSites;
    Rates2=MLELam(:,1);
    % %find the sites that intersect (should be all of them)
    [C,ia,ib]=intersect(Sites1,Sites2);
    
    K1 = log10(Rates1(ia)); % Log of the K values of the overlapped sites in data 1.
    K2 = log10(Rates2(ib)); % in data 41
    
    load(read_file2, 'AllDat', 'sites');
    [C2,ic,id]=intersect(sites, C);
    sz = 10;
    MAX_DEPTH1 = 20;
    MAX_DEPTH2 = 7;
    Read_Depth = sum(AllDat(ic, :, 1 : 2), 3);
    Read_Depth(Read_Depth(:, 1) > MAX_DEPTH1, 1) = MAX_DEPTH1;
    Read_Depth(Read_Depth(:, 2) > MAX_DEPTH2, 2) = MAX_DEPTH2;
    Read_Depth(Read_Depth(:, 3) > MAX_DEPTH2, 3) = MAX_DEPTH2;
    Read_Depth(Read_Depth(:, 4) > MAX_DEPTH2, 4) = MAX_DEPTH2;
    Read_Depth(:, 1) = Read_Depth(:, 1)./MAX_DEPTH1;
    Read_Depth(:, 2:4) = Read_Depth(:, 2:4)./MAX_DEPTH2;
    Read_Depth = uint8(round(Read_Depth * 10.0));
    N_Times = size(Read_Depth, 2);
    fig = figure(1);
    for ii = 1 : N_Times
        subplot(2, 2, ii);
        scatter(K1,K2, sz, Read_Depth(:, ii), 'filled');
        colorbar;
        xlabel(label1);
        ylabel(label2);
        xlim([-2, 1]);
        ylim([-2, 1]);
    end
    print(fig, figure_path, '-dpdf','-opengl','-r300');
    close;
end