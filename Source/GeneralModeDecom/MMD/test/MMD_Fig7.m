% This code shows an example when the number of components used for MMD
% (numGroupTest) is larger than the actual number of components (numGroup).
% The MMD is able to extract the numGroup components precisely and provide
% numGroupTest-numGroup components with O(eps) amplitude. The MMD is able
% to estimate the number of components with a proper thresholding parameter
% for the least energy of a component.
%
% Reference:
% "Multiresolution Mode Decomposition (MMD) for Adaptive Time Series Analysis"
% by Haizhao Yang.

if (1)
    close all;
    clear all;
    
    %% generate signal
    N = 2^14;      % sampling size
    NM = 0;       % noise level
    ns = NM*randn(1,N);
    
    x = [0:N-1]/N;
    t = x;
    amp = 0.006;%0.01
    F1 = 150;
    F2 = 220;
    
    sh1 = @(x) gen_shape2(x,3);
    sh2 = @(x) gen_shape2(x,2);
    
    numGroup = 2;
    numGroupTest = 4;
    
    xx = x + amp*sin(2*pi*x);
    insFreq(1,:) = (1+amp*2*pi*cos(2*pi*x))*F1;
    insAmp(1,:) = 1+0.1*sin(2*pi*xx)+0.2*cos(2*pi*xx);
    f1 = insAmp(1,:).*sh1(F1*xx);
    
    yy = x + amp*cos(2*pi*x);
    insFreq(2,:) = (1-amp*2*pi*sin(2*pi*x))*F2;
    insAmp(2,:) = 1+0.2*sin(2*pi*yy)+0.1*cos(2*pi*yy);
    f2 = insAmp(2,:).*sh2(F2*yy);
    
    insFreq(3,:) = 2*insFreq(1,:);
    insFreq(4,:) = 3*insFreq(2,:);
    
    insPhase(1,:) = (xx)*F1;
    insPhase(2,:) = (yy)*F2;
    insPhase(3,:) = 2*insPhase(1,:);
    insPhase(4,:) = 3*insPhase(2,:);
    
    N1 = 1; N2 = N; N3 = N/8;
    opt.maxiter = 200;
    opt.eps_error = 1e-6;
    opt.show = 0;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 2;
    opt.eps_diff = 1e-6;
    opt.ampErrBandWidth = 20;
    opt.numSweep = 10;
    opt.name = './results/MMD_fig7_reconcomp';
    opt.para.nknots = 20;
    opt.para.knotremoval_factor= 1.0001;
    opt.para.order = 3;
    opt.para.Ls = 1000;
    
    % test example: two components
    sig = f1 + f2  + ns;
    shapeTrue = cell(1,2);
    shapeTrue{1} = @(x) sh1(x);
    shapeTrue{2} = @(x) sh2(x);
    [shape,comp,Hcoef,flag,idx] = DeCom_MMD(sig,x,numGroupTest,ones(numGroupTest,N),insFreq,insPhase,opt);
    flag
    idx
    
    save('./results/MMD_fig7.mat','-v7.3');
    
    for i = 1:numGroup
        [trans_est_shape,min_error]=shape_phase_trans(shape{i}.s0,shapeTrue{i}(linspace(0,1,opt.para.Ls)));
        L = length(trans_est_shape);
        gd = 0:1/L:(1-1/L);
        pic = figure;plot(gd,trans_est_shape,'LineWidth',2);hold on;plot(gd,shapeTrue{i}(linspace(0,1,opt.para.Ls)),'LineWidth',2); hold off;
        legend('Est','True'); title([num2str(i),'th shape with error=',num2str(min_error)]);axis square;
        title('remove');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        tit = sprintf('./results/MMD_fig7_%d.fig',i);
        saveas(pic,tit);
        str = sprintf('./results/MMD_fig7_%d',i);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    end
    for i = numGroup+1:numGroupTest
        L = length(shape{i}.s0);
        gd = 0:1/L:(1-1/L);
        pic = figure;plot(gd,shape{i}.s0,'LineWidth',2);axis square;
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        tit = sprintf('./results/MMD_fig7_%d.fig',i);
        saveas(pic,tit);
        str = sprintf('./results/MMD_fig7_%d',i);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    end
    pic = figure;
    plot(x(N1:N3),comp{1}(N1:N3)); hold on;
    plot(x(N1:N3),f1(N1:N3)); hold off;
    xlabel('time');ylabel('signal intensity');legend('Est','True');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/MMD_fig7_checkcomp1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig7_checkcomp1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(x(N1:N3),comp{2}(N1:N3)); hold on;
    plot(x(N1:N3),f2(N1:N3)); hold off;
    xlabel('time');ylabel('signal intensity');legend('Est','True');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/MMD_fig7_checkcomp2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig7_checkcomp2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(x(N1:N3),comp{3}(N1:N3));
    xlabel('time');ylabel('signal intensity');legend('Est');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/MMD_fig7_checkcomp3.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig7_checkcomp3';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(x(N1:N3),comp{4}(N1:N3));
    xlabel('time');ylabel('signal intensity');legend('Est');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/MMD_fig7_checkcomp4.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig7_checkcomp4';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(x(N1:N3),f1(N1:N3)-comp{1}(N1:N3),'b');
    xlabel('time');ylabel('error');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/MMD_fig7_errcomp1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig7_errcomp1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(x(N1:N3),f2(N1:N3)-comp{2}(N1:N3),'b');
    xlabel('time');ylabel('error');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/MMD_fig7_errcomp2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig7_errcomp2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

