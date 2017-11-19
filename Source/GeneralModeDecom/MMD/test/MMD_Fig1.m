% This code generates Figure 1 of the paper
% "Multiresolution Mode Decomposition (MMD) for Adaptive Time Series Analysis"
% by Haizhao Yang.

if (1)
    close all;
    clear all;
    %% generate signal
    load 0009_8min.mat;
    N = 2^14;
    x = (0:N-1)/N;
    sig = signal.pleth.y(1:N)';
    
    numGroup = 2;
    opt.eps = 1e-3;
    opt.res = 0.25;
    opt.freq_range = [0 N/32];
    opt.NG = N;
    opt.dt = 1/N;
    opt.t_sc = 0.5;
    opt.NM = 0;
    opt.st = round([50 1200 ]/16/opt.res)*N/2^14;
    opt.ed = round([400 2000 ]/16/opt.res)*N/2^14;
    opt.num_select = numGroup;
    opt.red = 8;
    opt.C = 100;
    opt.rad = 1.5;
    opt.show = 0;
    [insFreq,insAmp,insPhase,comp_select] = insInfo(sig,opt);
    % find peaks
    peaks = zeros(size(comp_select));
    for cnt = 1:numGroup
        peaks(cnt,:) = peakDetection(comp_select(cnt,:),insFreq(cnt,:));
    end
    
    % correct phases
    insPhase = phaseShift(insPhase,peaks);
    for cnt = 1:2
        pic = figure;plot(x,insFreq(cnt,:),'LineWidth',2);
        axis square;
        title('remove');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        tit = sprintf('./results/fig1short_%d.fig',8+cnt-1);
        saveas(pic,tit);
        str = sprintf('./results/MMD_fig1short_%d',8+cnt-1);
        print(gcf, '-depsc', str);
        pic = figure;plot(x,insAmp(cnt,:),'LineWidth',2);
        axis square;
        title('remove');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        tit = sprintf('./results/fig1short_%d.fig',10+cnt-1);
        saveas(pic,tit);
        str = sprintf('./results/MMD_fig1short_%d',10+cnt-1);
        print(gcf, '-depsc', str);
    end
    save('./results/MMD_fig1.mat','-v7.3');
end

if (1)
    load ./results/MMD_fig1.mat;
    opt.maxiter = 200;
    opt.eps_error = 1e-6;
    opt.show = 0;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 2;
    opt.eps_diff = 1e-6;
    opt.ampErrBandWidth = 40;
    opt.numSweep = 10;
    
    opt.para.nknots = 20;
    opt.para.knotremoval_factor= 1.0001;
    opt.para.order = 3;
    opt.para.Ls = 1000;
    
    % test example: two components
    [shape,comp] = DeCom_MMD(sig,x,numGroup,ones(size(insAmp)),insFreq,insPhase,opt);
    
    for cnt1 = 1:numGroup
        shape{cnt1}.s0 = shape{cnt1}.s0/(norm(shape{cnt1}.s0)/sqrt(length(shape{cnt1}.s0)));
        LL = length(shape{cnt1}.s0);
        pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt1}.s0);
        title('remove');axis square;
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        tit = sprintf('./results/MMD_fig1short_%d.fig',i);
        saveas(pic,tit);
        str = sprintf('./results/MMD_fig1short_%d',i);
        print(gcf, '-depsc', str);
    end
    
    %% Modified to regenerate the signal to test if it is actually accurate enough
    Sig = comp{1} + comp{2};
    
    N1 = 1; N2 = N;
    Ntotal = length(signal.pleth.y);
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp{1}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/MMD_fig1_comp1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig1_comp1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp{2}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/MMD_fig1_comp2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig1_comp2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,Sig(N1:N2),'r');
    pbaspect([10 1 1]);set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    hold on
    plot((N1:N2)/Ntotal*480,sig(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    saveas(pic,'./results/MMD_fig1_sum.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig1_sum';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,sig(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/MMD_fig1_org.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig1_org';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,sig(N1:N2)-Sig(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/MMD_fig1_res.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig1_res';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end





