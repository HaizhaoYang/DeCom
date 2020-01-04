% This code generates Figure 1 to 3 of the paper
% "Multiresolution Mode Decomposition for Adaptive Time Series Analysis"
%
% by Haizhao Yang.

if (1)
    close all;
    clear all;
    %% generate signal
    fileName = {'0009_8min.mat','0015_8min.mat','0016_8min.mat','0018_8min.mat','0023_8min.mat','0028_8min.mat','0029_8min.mat',...
        '0030_8min.mat','0031_8min.mat','0032_8min.mat','0035_8min.mat','0038_8min.mat','0103_8min.mat','0104_8min.mat','0105_8min.mat',...
        '0115_8min.mat','0121_8min.mat','0122_8min.mat','0123_8min.mat','0125_8min.mat','0127_8min.mat','0128_8min.mat','0133_8min.mat',...
        '0134_8min.mat','0142_8min.mat','0147_8min.mat','0148_8min.mat','0149_8min.mat','0150_8min.mat','0309_8min.mat','0311_8min.mat',...
        '0312_8min.mat','0313_8min.mat','0322_8min.mat','0325_8min.mat','0328_8min.mat','0329_8min.mat','0330_8min.mat','0331_8min.mat',...
        '0332_8min.mat','0333_8min.mat','0370_8min.mat'};
    example = 3; % choose an example from the dataset above
    filePath = sprintf('./Applications/MedicalSci/Data/TBME2013-PPGRR-Benchmark_R3/data/%s',fileName{example});
    load(filePath);
    N = 2^14;
    x = (0:N-1)/N;
    sig = signal.pleth.y(N+1:2*N)';
    
    numGroup = 2;
    opt.eps = 1e-3;
    opt.res = 0.25;
    opt.freq_range = [0 N/32];
    [posl,posh,shat] = freqInfoPPG(sig);
    posl = posl/opt.res;
    posh = posh/opt.res;
    opt.st = round([posl*0.75 posh*0.75 ]);
    opt.ed = round([posl*1.25 posh*1.25 ]);
    opt.NG = N;
    opt.dt = 1/N;
    opt.t_sc = 0.5;
    opt.NM = 0;
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
end

N1 = 1; N2 = N;
Ntotal = length(signal.pleth.y);
pic = figure;
plot((N1:N2)/Ntotal*480,sig(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
title('Original PPG signal');
saveas(pic,'./results/MMD_fig1_org.fig');
set(gca, 'FontSize', 16);
b=get(gca);
set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
str = './results/MMD_fig1_org';
print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);

% perform MMD

if (1)
    opt.maxiter = 200;
    opt.eps_error = 1e-6;
    opt.show = 0;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 2;
    opt.eps_diff = 1e-6;
    opt.ampErrBandWidth = 40;
    opt.numSweep = 10;
    
    switch opt.shapeMethod
        case 1
            opt.para.Ls=1000;
            opt.para.bandWidth = 10;
            opt.para.diffeoMethod = 'nufft';
        case 2
            opt.para.nknots = 10;
            opt.para.knotremoval_factor= 1.0001;
            opt.para.order = 3;
            opt.para.Ls = 1000;
    end
    
    % test example: two components
    [shape,comp] = DeCom_MMD(sig,x,numGroup,ones(size(insAmp)),insFreq,insPhase,opt);
    
    %% Modified to regenerate the signal to test if it is actually accurate enough
    Sig = comp{1} + comp{2};
    opt.show = 1;
    insInfo(sig-Sig,opt);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp{1}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    title('Component 1');
    saveas(pic,'./results/MMD_fig1_comp1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig1_comp1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp{2}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    title('Component 2');
    saveas(pic,'./results/MMD_fig1_comp2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig1_comp2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,sig(N1:N2)-Sig(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    title('Residual');
    saveas(pic,'./results/MMD_fig1_res.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig1_res';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

% perform GMD

if (1)
    opt.maxiter = 200;
    opt.eps_error = 1e-6;
    opt.show = 0;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 2;
    opt.eps_diff = 1e-6;
    opt.ampErrBandWidth = 0;
    opt.numSweep = 1;
    
    switch opt.shapeMethod
        case 1
            opt.para.Ls=1000;
            opt.para.bandWidth = 10;
            opt.para.diffeoMethod = 'nufft';
        case 2
            opt.para.nknots = 10;
            opt.para.knotremoval_factor= 1.0001;
            opt.para.order = 3;
            opt.para.Ls = 1000;
    end
    
    % test example: two components
    [shape1,comp1] = DeCom_MMD(sig,x,numGroup,insAmp,insFreq,insPhase,opt);
    
    %% Modified to regenerate the signal to test if it is actually accurate enough
    Sig1 = comp1{1} + comp1{2};
    save('./results/MMD_fig1.mat','-v7.3');
    
    insInfo(sig-Sig1,opt);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp1{1}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    title('Component 1');
    saveas(pic,'./results/MMD_fig1_comp1_p2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig1_comp1_p2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,comp1{2}(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    title('Component 2');
    saveas(pic,'./results/MMD_fig1_comp2_p2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig1_comp2_p2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot((N1:N2)/Ntotal*480,sig(N1:N2)-Sig1(N1:N2),'b'); axis tight; xlabel('time');ylabel('signal intensity');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);%axis([N1/Ntotal*480, N2/Ntotal*480,-15,15]);
    title('Residual');
    saveas(pic,'./results/MMD_fig1_res_p2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig1_res_p2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end
