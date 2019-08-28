% This code generates Figure 11-15 of the paper
%
% "A Fast Algorithm for Multiresolution Mode Decomposition"
% by Haizhao Yang.

if (1)
    close all;
    clear all;
    
    %% generate signal
    N = 2^16;  %18    %% sampling size
    NM = 0;       %% noise level
    ns = NM*randn(1,N);
    
    x = [0:N-1]/N;
    t = x;
    amp = 0.006;%0.01
    F1 = 150;
    F2 = 220;
    
    sh1 = @(x) gen_shape2(x,3);
    sh2 = @(x) gen_shape2(x,2);
    
    numGroup = 2;
    
    xx = x + amp*sin(2*pi*x);
    insFreq(1,:) = (1+amp*2*pi*cos(2*pi*x))*F1;
    a1sin = 0.1*sin(2*pi*xx);
    a1cos = 0.2*cos(2*pi*xx); % error
    insAmp(1,:) = 1+a1sin+a1cos;
    f1 = insAmp(1,:).*sh1(F1*xx); %%%%%ECG gen_shape(FF*xx,2)
    
    yy = x + amp*cos(2*pi*x);
    insFreq(2,:) = (1-amp*2*pi*sin(2*pi*x))*F2;
    a2sin = 0.2*sin(2*pi*yy); %error
    a2cos = 0.1*cos(2*pi*yy);
    insAmp(2,:) = 1+a2sin+a2cos;
    f2 = insAmp(2,:).*sh2(F2*yy); %%%%%%
    
    insPhase(1,:) = (xx)*F1;
    insPhase(2,:) = (yy)*F2;
    
    N1 = 1; N2 = N; N3 = N/8;
    if (0)
        % show estimation target
        pic = figure;
        plot(x(N1:N2),sh1(F1*xx(N1:N2)),'b'); axis tight; xlabel('time');ylabel('signal intensity');
        pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);
        saveas(pic,'./results/RDSA_fig11_2_truecomp1_1.fig');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        str = './results/RDSA_fig11_2_truecomp1_1';
        print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
        
        
        pic = figure;
        plot(x(N1:N2),sh2(F2*yy(N1:N2)),'b'); axis tight; xlabel('time');ylabel('signal intensity');
        pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);
        saveas(pic,'./results/RDSA_fig11_2_truecomp2_1.fig');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        str = './results/RDSA_fig11_2_truecomp2_1';
        print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
        
        pic = figure;
        plot(x(N1:N2),a1sin(N1:N2).*sh1(F1*xx(N1:N2)),'b'); axis tight; xlabel('time');ylabel('signal intensity');
        pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);
        saveas(pic,'./results/RDSA_fig11_2_truecompsin1_2.fig');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        str = './results/RDSA_fig11_2_truecompsin1_1';
        print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
        
        
        pic = figure;
        plot(x(N1:N2),a2cos(N1:N2).*sh2(F2*yy(N1:N2)),'b'); axis tight; xlabel('time');ylabel('signal intensity');
        pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);
        saveas(pic,'./results/RDSA_fig11_2_truecompcos2_2.fig');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        str = './results/RDSA_fig11_2_truecompcos2_1';
        print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
        
        
        pic = figure;
        plot(x(N1:N2),a1cos(N1:N2).*sh1(F1*xx(N1:N2)),'b'); axis tight; xlabel('time');ylabel('signal intensity');
        pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);
        saveas(pic,'./results/RDSA_fig11_2_truecompcos1_2.fig');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        str = './results/RDSA_fig11_2_truecompcos1_1';
        print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
        
        
        pic = figure;
        plot(x(N1:N2),a2sin(N1:N2).*sh2(F2*yy(N1:N2)),'b'); axis tight; xlabel('time');ylabel('signal intensity');
        pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]);
        saveas(pic,'./results/RDSA_fig11_2_truecompsin2_2.fig');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        str = './results/RDSA_fig11_2_truecompsin2_1';
        print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    end
    
    % test example: two components
    numGroup = 2;
    sig = f1 + f2  + ns;
    shapeTrue = cell(1,2);
    shapeTrue{1} = @(x) sh1(x);
    shapeTrue{2} = @(x) sh2(x);
    
    if (1)
        % estimate instantaneous frequency and phase
        opt.eps = 1e-3;
        opt.res = 0.05;
        opt.freq_range = [0 N/128];
        opt.NG = N/16;
        opt.dt = 1/N;
        opt.t_sc = 0.75;
        opt.NM = 0;
        opt.st = round([500 800 ]*5);
        opt.ed = round([700 1000 ]*5);
        opt.num_select = numGroup;
        opt.red = 8;
        opt.C = 100;
        opt.rad = 1;
        opt.show = 1;
        [insFreq_est,insAmp_est,insPhase_est,comp_est] = insInfo(sig,opt);
        insFreq_est = insFreq_est*2; insPhase_est = insPhase_est*2;
        % find peaks
        peaks = zeros(size(comp_est));
        for cnt = 1:numGroup
            peaks(cnt,:) = peakDetection(comp_est(cnt,:),insFreq_est(cnt,:));
        end
        % correct phases
        insPhase_est = phaseShift(insPhase_est,peaks);
        
        % plot instantaneous property estimation
        
        pic = figure;axis square;
        plot((insPhase(1,:)-insPhase_est(1,:)-insPhase(1,1)),'b');
        hold on; plot((insPhase(2,:)-insPhase_est(2,:)-insPhase(2,1)),'r');
        legend('error 1','error 2');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        str = sprintf('./results/RDSA_fig11_2_IP'); saveas(pic,str);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
        pic = figure;axis square;
        h(1) = plot(x,insFreq_est(1,:),'b','LineWidth',2);
        hold on;
        h(2) = plot(x,insFreq(1,:),'r','LineWidth',2);
        hold off;
        legend('Est','True');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        str = sprintf('./results/RDSA_fig11_2_IF1'); saveas(pic,str);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
        pic = figure;axis square;
        h(1) = plot(x,insFreq_est(2,:),'b','LineWidth',2);
        hold on;
        h(2) = plot(x,insFreq(2,:),'r','LineWidth',2);
        hold off;
        legend('Est','True');
        set(gca, 'FontSize', 16);
        b=get(gca);
        set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
        str = sprintf('./results/RDSA_fig11_2_IF2'); saveas(pic,str);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
    end
    
    opt.maxiter = 200;
    opt.eps_error = 1e-6;
    opt.show = 0;
    opt.iterStyle = 'GS';
    opt.shapeMethod = 1;
    opt.eps_diff = 1e-6;
    opt.ampErrBandWidth = 20;
    opt.name = './results/RDSA_fig11_2_reconcomp';
    opt.numSweep = 10;
    switch opt.shapeMethod
        case 1
            opt.para.Ls=2000;
            opt.para.bandWidth = 2000;
        case 2
            opt.para.nknots = 20;
            opt.para.knotremoval_factor= 1.0001;
            opt.para.order = 3;
            opt.para.Ls = 1000;
    end
    
    % estimate shapes
    insPhase_est(1,:) = insPhase_est(1,:) + insPhase(1,1);% just for the purpose of better visualization
    insPhase_est(2,:) = insPhase_est(2,:) + insPhase(2,1);
    [shape,comp,Hcoef] = DeCom_MMD(sig,x,numGroup,ones(size(insAmp)),insFreq_est,insPhase_est,opt);
    
    save('./results/RDSA_fig11_2.mat','-v7.3');
end

if (1)
    load ./results/RDSA_fig11_2.mat;
    
    close all;
    
    Tcoef = cell(1,2);
    Tcoef{1}.a0 = 1; Tcoef{2}.a0 = 1;
    Tcoef{1}.a = zeros(1,10); Tcoef{1}.a(1) = 0.2;
    Tcoef{1}.b = zeros(1,10); Tcoef{1}.b(1) = 0.1;
    Tcoef{2}.a = zeros(1,10); Tcoef{2}.a(1) = 0.1;
    Tcoef{2}.b = zeros(1,10); Tcoef{2}.b(1) = 0.2;
    
    for cnt = 1:numGroup
        LL = length(shape{cnt}.s0);
        pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt}.s0*Hcoef{cnt}.a0);
        hold on; plot(0:1/LL:(1-1/LL),Tcoef{cnt}.a0*shapeTrue{cnt}(linspace(0,1,LL))/sqrt(sum(abs(shapeTrue{cnt}(linspace(0,1,LL))).^2)/LL),'LineWidth',2); hold off;
        if cnt == 1, legend('Est','True','Location','Southeast'); else legend('Est','True','Location','Northeast'); end;
        title('remove');axis square;
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/RDSA_fig11_2_%d_shape_%d.fig',cnt,0);
        saveas(pic,tit);
        str = sprintf('./results/RDSA_fig11_2_%d_shape_%d',cnt,0);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
        LL = length(shape{cnt}.scn{1});
        pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt}.scn{1}*Hcoef{cnt}.a(1));
        hold on; plot(0:1/LL:(1-1/LL),Tcoef{cnt}.a(1)*shapeTrue{cnt}(linspace(0,1,LL))/sqrt(sum(abs(shapeTrue{cnt}(linspace(0,1,LL))).^2)/LL),'LineWidth',2); hold off;
        if cnt == 1, legend('Est','True','Location','Southeast'); else legend('Est','True','Location','Northeast'); end;
        title('remove');axis square;
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/RDSA_fig11_2_%d_shapeCn_%d.fig',cnt,1);
        saveas(pic,tit);
        str = sprintf('./results/RDSA_fig11_2_%d_shapeCn_%d',cnt,1);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
        LL = length(shape{cnt}.ssn{1});
        pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt}.ssn{1}*Hcoef{cnt}.b(1));
        hold on; plot(0:1/LL:(1-1/LL),Tcoef{cnt}.b(1)*shapeTrue{cnt}(linspace(0,1,LL))/sqrt(sum(abs(shapeTrue{cnt}(linspace(0,1,LL))).^2)/LL),'LineWidth',2); hold off;
        if cnt == 1, legend('Est','True','Location','Southeast'); else legend('Est','True','Location','Northeast'); end;
        title('remove');axis square;
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/RDSA_fig11_2_%d_shapeSn_%d.fig',cnt,1);
        saveas(pic,tit);
        str = sprintf('./results/RDSA_fig11_2_%d_shapeSn_%d',cnt,1);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
        LL = length(shape{cnt}.scn{2});
        pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt}.scn{2}*Hcoef{cnt}.a(2));
        hold on; plot(0:1/LL:(1-1/LL),Tcoef{cnt}.a(2)*shapeTrue{cnt}(linspace(0,1,LL))/sqrt(sum(abs(shapeTrue{cnt}(linspace(0,1,LL))).^2)/LL),'LineWidth',2); hold off;
        if cnt == 1, legend('Est','True','Location','Southeast'); else legend('Est','True','Location','Northeast'); end;
        title('remove');axis square;
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/RDSA_fig11_2_%d_shapeCn_%d.fig',cnt,-1);
        saveas(pic,tit);
        str = sprintf('./results/RDSA_fig11_2_%d_shapeCn_%d',cnt,-1);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
        LL = length(shape{cnt}.ssn{2});
        pic = figure;plot(0:1/LL:(1-1/LL),shape{cnt}.ssn{2}*Hcoef{cnt}.b(2));
        hold on; plot(0:1/LL:(1-1/LL),Tcoef{cnt}.b(2)*shapeTrue{cnt}(linspace(0,1,LL))/sqrt(sum(abs(shapeTrue{cnt}(linspace(0,1,LL))).^2)/LL),'LineWidth',2); hold off;
        if cnt == 1, legend('Est','True','Location','Southeast'); else legend('Est','True','Location','Northeast'); end;
        title('remove');axis square;
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/RDSA_fig11_2_%d_shapeSn_%d.fig',cnt,-1);
        saveas(pic,tit);
        str = sprintf('./results/RDSA_fig11_2_%d_shapeSn_%d',cnt,-1);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
        
    end
    
    
    for i = 1:numGroup
        [trans_est_shape,min_error]=shape_phase_trans(shape{i}.s0,shapeTrue{i}(linspace(0,1,opt.para.Ls)));
        L = length(trans_est_shape);
        gd = 0:1/L:(1-1/L);
        pic = figure;plot(gd,shapeTrue{i}(linspace(0,1,opt.para.Ls)),'LineWidth',2);
        title([num2str(i),'th shape with error=',num2str(min_error)]);axis square;
        title('remove');
        set(gca, 'FontSize', 32);
        b=get(gca);
        set(b.XLabel, 'FontSize', 32);set(b.YLabel, 'FontSize', 32);set(b.ZLabel, 'FontSize', 32);set(b.Title, 'FontSize', 32);
        tit = sprintf('./results/RDSA_fig11_2_%d.fig',i);
        saveas(pic,tit);
        str = sprintf('./results/RDSA_fig11_2_%d',i);
        print(gcf, '-depsc', str);   command = sprintf('epstopdf %s.eps',str);      system(command);
    end
    pic = figure;
    plot(x(N1:N3),comp{1}(N1:N3)); hold on;
    plot(x(N1:N3),f1(N1:N3)); hold off;
    xlabel('time');ylabel('signal intensity');legend('Est','True');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/RDSA_fig11_2_checkcomp1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig11_2_checkcomp1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(x(N1:N3),comp{2}(N1:N3)); hold on;
    plot(x(N1:N3),f2(N1:N3)); hold off;
    xlabel('time');ylabel('signal intensity');legend('Est','True');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/RDSA_fig11_2_checkcomp2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig11_2_checkcomp2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(x(N1:N3),f1(N1:N3)-comp{1}(N1:N3),'b');
    xlabel('time');ylabel('error');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/RDSA_fig11_2_errcomp1.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig11_2_errcomp1';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    pic = figure;
    plot(x(N1:N3),f2(N1:N3)-comp{2}(N1:N3),'b');
    xlabel('time');ylabel('error');
    pbaspect([10 1 1]); set(pic, 'Position', [200, 200, 1200, 200]); axis tight;
    saveas(pic,'./results/RDSA_fig11_2_errcomp2.fig');
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    str = './results/RDSA_fig11_2_errcomp2';
    print(gcf, '-depsc', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

