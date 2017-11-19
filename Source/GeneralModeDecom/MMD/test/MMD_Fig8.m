% This code shows an example when the number of components used for MMD
% (numGroupTest) is smaller than the actual number of components (numGroup).
% The MMD is able to extract the numGroupTest components approximately.
% The MMD is able to warn that the actual number of components should be
% larger.
%
% Reference:
% "Multiresolution Mode Decomposition (MMD) for Adaptive Time Series Analysis"
% by Haizhao Yang.

if (1)
    close all;
    clear all;
    
    %% generate signal
    change_amplitude = 1;
    N = 2^16;%18;      %% sampling size
    dif = 0.05;
    
    x = [0:N-1]/N;
    t = x;
    amp = 0.01;%0.01
    FF = 160;%160
    
    numGroup = 4;
    numGroupTest = 2;
    insFreq = zeros(numGroup,N);
    insAmp = zeros(numGroup,N);
    insPhase = zeros(numGroup,N);
    f = cell(1,numGroup);
    
    xx = x + amp*sin(2*pi*x);
    insFreq(1,:) = (1+amp*2*pi*cos(2*pi*x))*FF;
    f{1} = zeros(1,N);
    insAmp(1,:) = 1+0.05*sin(4*pi*xx);
    f{1} = insAmp(1,:).*gen_shape2(FF*xx,3); %%%%%ECG gen_shape(FF*xx,2)
    
    yy = x + dif + amp*sin(2*pi*(x+dif));
    insFreq(2,:) = (1+amp*2*pi*cos(2*pi*(x+dif)))*FF;
    f{2} = zeros(1,N);
    insAmp(2,:) = 1+0.05*sin(2*pi*yy);
    f{2} = insAmp(2,:).*gen_shape(FF*yy,2); %%%%%%
    
    
    zz = x + dif*2 + amp*sin(2*pi*(x+dif*2));
    insFreq(3,:) = (1+amp*2*pi*cos(2*pi*x+dif*2))*FF;
    f{3} = zeros(1,N);
    insAmp(3,:) = 1+0.05*sin(2*pi*zz);
    f{3} = insAmp(3,:).*gen_shape(FF*zz,1);%%%%%%%%%%%%
    
    ww = x + dif*3 + amp*sin(2*pi*(x+dif*3));
    insFreq(4,:) = (1+amp*2*pi*cos(2*pi*x+dif*3))*FF;
    f{4} = zeros(1,N);
    insAmp(4,:) = 1+0.05*sin(2*pi*ww);
    f{4} = insAmp(4,:).*gen_shape(FF*ww,3);%%%%%%%%%%
    
    insPhase(1,:) = (xx)*FF;
    insPhase(2,:) = (yy)*FF;
    insPhase(3,:) = (zz)*FF;
    insPhase(4,:) = (ww)*FF;
    
    
    opt.maxiter = 200;
    opt.eps_error = 1e-8;
    opt.show = 0;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 2;
    opt.eps_diff = 1e-6;
    opt.ampErrBandWidth = 20;
    opt.lowestCompEng = 1e-1;
    opt.decayBndExpCoef = 0.01;
    opt.numSweep = 10;
    opt.para.nknots = 20;
    opt.para.knotremoval_factor= 1.0001;
    opt.para.order = 3;
    opt.para.Ls = 1000;
    
    sig = f{1} + f{2} +f{3} + f{4};
    shapeTrue = cell(1,4);
    shapeTrue{1} = @(x) gen_shape2(x,3);
    shapeTrue{2} = @(x) gen_shape(x,2);
    shapeTrue{3} = @(x) gen_shape(x,1);
    shapeTrue{4} = @(x) gen_shape(x,3);
    [shape,comp,Hcoef,flag,idx] = DeCom_MMD(sig,x,numGroupTest,ones(numGroupTest,N),insFreq(1:numGroupTest,:),insPhase(1:numGroupTest,:),opt);
    flag
    idx
    
    save('./results/MMD_fig8.mat','-v7.3');
    
    N1 = 1; N2 = N; N3 = N/8;
    for i = 1:numGroupTest
        [trans_est_shape,min_error]=shape_phase_trans(shape{i}.s0,shapeTrue{i}(linspace(0,1,opt.para.Ls)));
        L = length(trans_est_shape);
        gd = 0:1/L:(1-1/L);
        pic = figure;plot(gd,trans_est_shape,'LineWidth',2);hold on;plot(gd,shapeTrue{i}(linspace(0,1,opt.para.Ls)),'LineWidth',2); hold off;
        legend('Est','True'); title([num2str(i),'th shape with error=',num2str(min_error)]);axis square;
        title('remove');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        tit = sprintf('./results/MMD_fig8_%d.fig',i);
        saveas(pic,tit);
        str = sprintf('./results/MMD_fig8_%d',i);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    end
    pic = figure;
    plot(x(N1:N3),comp{1}(N1:N3)); hold on;
    plot(x(N1:N3),f{1}(N1:N3)); hold off;
    xlabel('time');ylabel('signal intensity');legend('Est','True');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/MMD_fig8_checkcomp1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig8_checkcomp1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(x(N1:N3),comp{2}(N1:N3)); hold on;
    plot(x(N1:N3),f{2}(N1:N3)); hold off;
    xlabel('time');ylabel('signal intensity');legend('Est','True');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/MMD_fig8_checkcomp2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig8_checkcomp2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(x(N1:N3),f{1}(N1:N3)-comp{1}(N1:N3),'b');
    xlabel('time');ylabel('error');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/MMD_fig8_errcomp1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig8_errcomp1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(x(N1:N3),f{2}(N1:N3)-comp{2}(N1:N3),'b');
    xlabel('time');ylabel('error');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/MMD_fig8_errcomp2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/MMD_fig8_errcomp2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
end


