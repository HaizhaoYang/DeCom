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
    
    opt.maxiter = 200;
    opt.eps_error = 1e-6;
    opt.show = 0;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 2;
    opt.eps_diff = 1e-6;
    opt.ampErrBandWidth = 30;
    opt.numSweep = 50;
    
    opt.para.nknots = 20;
    opt.para.knotremoval_factor= 1.0001;
    opt.para.order = 3;
    opt.para.Ls = 1000;
    
    % test example: two components
    [shape,comp] = DeCom_MMD(sig,x,numGroup,ones(size(insAmp)),insFreq,insPhase,opt);
    save('./results/MMD_fig0.mat','-v7.3');
end

if (1)
    load ./results/MMD_fig0.mat;
    close all;
    opt.eps = 1e-1;
    N2 = 3*N/4;    N1 = 1+N/2;
    Ntotal = length(signal.pleth.y);
    opt.res = 0.2;
    opt.t_sc = 0.8;
    T_f = ss_wp1_fwd(comp{2},1,1,1,x,opt.NG,opt.freq_range(2)*2,opt.freq_range(1),opt.rad,1,opt.t_sc,opt.red,opt.eps,opt.res);
    st = 750;
    ed = 1000;
    pic = figure;imagesc([N1, N2]/Ntotal*480,[st ed]*opt.res/(N/Ntotal*480),real(T_f(st:ed,N1:N2)));
    xlabel('time (Second)');ylabel('frequency (Hz)');axis xy;
    % pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);
    colormap (1-gray);
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = sprintf('./results/MMD_fig0_1');
    tit = [str,'.fig'];
    saveas(pic,tit);
    print(gcf, '-depsc', str);    command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    st = 400;
    ed = 550;
    pic = figure;imagesc([N1, N2]/Ntotal*480,[st ed]*opt.res/(N/Ntotal*480),real(T_f(st:ed,N1:N2)));
    xlabel('time (Second)');ylabel('frequency (Hz)');axis xy;
    % pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);
    colormap (1-gray);
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = sprintf('./results/MMD_fig0_1');
    tit = [str,'.fig'];
    saveas(pic,tit);
    print(gcf, '-depsc', str);    command = sprintf('epstopdf %s.eps',str);      system(command);
    
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp{2}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/MMD_fig0_comp2_2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig0_comp2_2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    opt.ampErrBandWidth = 0;
    opt.maxiter = 2;
    [shape,comp] = DeCom_MMD(sig,x,numGroup,ones(size(insAmp)),insFreq,insPhase,opt);
    
    T_f = ss_wp1_fwd(comp{2},1,1,1,x,opt.NG,opt.freq_range(2)*2,opt.freq_range(1),opt.rad,1,opt.t_sc,opt.red,opt.eps,opt.res);
    st = 750;
    ed = 1000;
    pic = figure;imagesc([N1, N2]/Ntotal*480,[st ed]*opt.res/(N/Ntotal*480),real(T_f(st:ed,N1:N2)));
    xlabel('time (Second)');ylabel('frequency (Hz)');axis xy;
    %pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);
    colormap (1-gray);
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = sprintf('./results/MMD_fig0_2');
    tit = [str,'.fig'];
    saveas(pic,tit);
    print(gcf, '-depsc', str);    command = sprintf('epstopdf %s.eps',str);      system(command);
    
    st = 400;
    ed = 550;
    pic = figure;imagesc([N1, N2]/Ntotal*480,[st ed]*opt.res/(N/Ntotal*480),real(T_f(st:ed,N1:N2)));
    xlabel('time (Second)');ylabel('frequency (Hz)');axis xy;
    %pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);
    colormap (1-gray);
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = sprintf('./results/MMD_fig0_2');
    tit = [str,'.fig'];
    saveas(pic,tit);
    print(gcf, '-depsc', str);    command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp{2}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    saveas(pic,'./results/MMD_fig0_comp2_3.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig0_comp3';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end





