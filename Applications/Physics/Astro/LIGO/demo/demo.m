clear all;
close all;

%% data information
% H-H1_LOSC_4_V1-1126257414-4096.hdf5  4096 seconds,  4096 Hz, Hanford
% H-H1_LOSC_4_V1-1126259446-32.hdf5    32   seconds,  4096 Hz, Hanford
% H-H1_LOSC_16_V1-1126259446-32.hdf5   32   seconds, 16384 Hz, Hanford
% H-H1_LOSC_16_V1-1126257414-4096.hdf5 4096 seconds, 16384 Hz, Hanford
% L-L1_LOSC_4_V1-1126257414-4096.hdf5  4096 seconds,  4096 Hz, Livingston
% L-L1_LOSC_4_V1-1126259446-32.hdf5    32   seconds,  4096 Hz, Livingston
% L-L1_LOSC_16_V1-1126259446-32.hdf5   32   seconds, 16384 Hz, Livingston
% L-L1_LOSC_16_V1-1126257414-4096.hdf5 4096 seconds, 16384 Hz, Livingston

%% set up data
% choose data set
dSet = 1;
switch dSet
    case 1
        FILENAME = 'H-H1_LOSC_4_V1-1126257414-4096.hdf5';
        sampleRate = 4096; timeSt = 0; timeEnd = 4096;
    case 2
        FILENAME = 'H-H1_LOSC_4_V1-1126259446-32.hdf5';
        sampleRate = 4096; timeSt = 0; timeEnd = 32;
    case 3
        FILENAME = 'H-H1_LOSC_16_V1-1126257414-4096.hdf5';
        sampleRate = 16384; timeSt = 0; timeEnd = 4096;
    case 4
        FILENAME = 'H-H1_LOSC_16_V1-1126259446-32.hdf5';
        sampleRate = 16384; timeSt = 0; timeEnd = 32;
    case 5
        FILENAME = 'L-L1_LOSC_4_V1-1126257414-4096.hdf5';
        sampleRate = 4096; timeSt = 0; timeEnd = 4096;
    case 6
        FILENAME = 'L-L1_LOSC_4_V1-1126259446-32.hdf5';
        sampleRate = 4096; timeSt = 0; timeEnd = 32;
    case 7
        FILENAME = 'L-L1_LOSC_16_V1-1126257414-4096.hdf5';
        sampleRate = 16384; timeSt = 0; timeEnd = 4096;
    case 8
        FILENAME = 'L-L1_LOSC_16_V1-1126259446-32.hdf5';
        sampleRate = 16384; timeSt = 0; timeEnd = 32;
end
% read data set
FILEINFO = hdf5info(FILENAME)
data = hdf5read(FILEINFO.GroupHierarchy.Groups(3).Datasets(1));
figure;plot(data);
datahat = fftshift(fft(data));
N0 = length(data);
figure;plot(-N0/2:(N0/2-1),abs(datahat));

N2 = 2^floor(log2(length(data)))
display('Length of the data segment:');
N = min(2^14,2^floor(log2(length(data))))

% pre-processing
x = [0:N-1]/N;
if 0
    sig = data(end-N2+1:end);
    ds = N2/N;
    sig = sig(1:ds:end); % subsample the data
    sampleRate = sampleRate/ds; % update the sampling rate
else
    sig = data(round(N0/2)-N/2+1:round(N0/2)+N/2);
end
sig = sig(:).';
sig = sig/max(sig);


% assume the signal is defined on the time interval [timeSt,timeEnd]
timeEnd = min(N/sampleRate,timeEnd);
% assume that we want to see the time-frequency plane in the frequency
% range of [lowBound,upBound]
upBound =  1100; % Hz
lowBound = 50;  % Hz


is_real = 1; % the signal is real

%% set up parameters
res = 1;
NG = N/2;
is_unif = 1;
typeNUFFT = 1;
is_cos = 1;
epsl = 1e-1;
red = 50;

freq_range = [lowBound,upBound];% assuming the data is defined on the time domain [0,1]
t_sc = 0.75;
rad = 1;

%% apply SS wave packet transform
tL = (timeEnd-timeSt)/4;
tC = (timeEnd+timeSt)/2;
[T_f coef kk] = ss_wp1_fwd(sig,is_real,is_unif,typeNUFFT,x,NG,upBound,lowBound,rad,is_cos,t_sc,red,1e-8,res,0);
pic = figure;imagesc([tC-tL tC+tL],[lowBound upBound]/res,(T_f(:,end/4+1:3*end/4)));colorbar off;title('SST');axis square;axis xy;
xlabel('x'); ylabel('k');
pic = figure;imagesc([tC-tL tC+tL],[lowBound upBound]/res,log2(T_f(:,end/4+1:3*end/4)));colorbar off;title('SST');axis square;axis xy;
xlabel('x'); ylabel('k');

%% second-order synchrosqueezing
N = min(2^13,2^floor(log2(length(data))))

% pre-processing
x = [0:N-1]/N;
if 0
    sig = data(end-N2+1:end);
    ds = N2/N;
    sig = sig(1:ds:end); % subsample the data
    sampleRate = sampleRate/ds; % update the sampling rate
else
    sig = data(round(N0/2)-N/2+1:round(N0/2)+N/2);
end
sig = sig(:).';
sig = sig/max(sig);


% assume the signal is defined on the time interval [timeSt,timeEnd]
timeEnd = min(N/sampleRate,timeEnd);
% Parameters
gamma = 0.0000001;
Nfft = N;
[~,~,~,VSST1] = sst2_new(sig,0.2,Nfft,gamma);
figure;imagesc([tC-tL tC+tL],[lowBound upBound],abs(VSST1(lowBound:upBound,:)));title('new VSST');colorbar off;title('SST2');axis square;axis xy;
xlabel('x'); ylabel('k');
figure;imagesc([tC-tL tC+tL],[lowBound upBound],log2(abs(VSST1(lowBound:upBound,:))));title('new VSST');colorbar off;title('SST2');axis square;axis xy;
xlabel('x'); ylabel('k');

%% only squeeze the largest coefficient
if 0
    ccc = wp1_fwd(x,is_real,is_unif,typeNUFFT,x,NG,upBound,lowBound,rad,is_cos, t_sc, red, 1);
    aaa = wp1_ext(x,is_real,is_unif,typeNUFFT,x,NG,upBound,lowBound,rad,is_cos, t_sc, red, 1);
    
    if is_real
        fqscale = [lowBound upBound];
    else
        fqscale = [-upBound upBound];
    end
    num_grid = ceil((fqscale(2)-fqscale(1))/res);
    res = (fqscale(2)-fqscale(1))/(num_grid-1);
    %grid = fqscale(1):res:fqscale(2);
    
    EXT = 10^10;
    ncl = cell(1,red);
    nclp = cell(1,red);
    NG = length(ccc{1}{1}{1});
    EXT = 10^10;
    kk = cell(1,red);
    T_f = zeros(num_grid,NG);
    coef = cell(1,red);
    loc = zeros(1,NG);
    val = zeros(1,NG);
    for cntred = 1:red
        ncl{cntred} = numel(ccc{cntred});
        nclp{cntred} = 0;
        for cnt = 1:ncl{cntred}
            nclp{cntred} = nclp{cntred} + numel(ccc{cntred}{cnt});
        end
        
        if is_real
            cc = zeros(nclp{cntred}/2,NG);
            aa = zeros(nclp{cntred}/2,NG);
            cnt_nclp = 1;
            for cnt = 1:ncl{cntred}
                cc(cnt_nclp,:) = ccc{cntred}{cnt}{1};
                aa(cnt_nclp,:) = aaa{cntred}{cnt}{1};
                cnt_nclp = cnt_nclp + 1;
            end
        else
            cc = zeros(nclp{cntred}-1,NG);
            aa = zeros(nclp{cntred}-1,NG);
            cnt_nclp1 = 1;
            cnt_nclp2 = 1;
            for cnt = 1:ncl{cntred}
                for cnt2 = 1:2
                    if cnt2 == 1
                        cc(cnt_nclp1+nclp{cntred}/2-1,:) = ccc{cntred}{cnt}{cnt2};
                        aa(cnt_nclp1+nclp{cntred}/2-1,:) = aaa{cntred}{cnt}{cnt2};
                        cnt_nclp1 = cnt_nclp1 + 1;
                    else if cnt ~= 1
                            cc(-cnt_nclp2+nclp{cntred}/2,:) = ccc{cntred}{cnt}{cnt2};
                            aa(-cnt_nclp2+nclp{cntred}/2,:) = aaa{cntred}{cnt}{cnt2};
                            cnt_nclp2 = cnt_nclp2 + 1;
                        end
                    end
                end
            end
        end
        
        szc = size(cc);
        gud = find(abs(cc)>epsl);
        aa=real(aa(gud)./(2*pi*i*cc(gud)));
        if is_real
            good = find(aa<=upBound & aa>=lowBound);
        else
            good = find(aa<=upBound & aa>=-upBound);
        end
        
        kk{cntred} = repmat(EXT,szc);
        kk{cntred}(gud(good)) = aa(good);
        
        coef{cntred} = zeros(size(cc));
        coef{cntred}(gud(good)) = cc(gud(good));
    end
    for cntred = 1:red
        temp = [val;abs(coef{cntred})];
        [maxval,maxloc] = max(temp);
        tempkk = [loc;kk{cntred}];
        val = maxval;
        for cnt = 1:length(loc)
            loc(cnt) = tempkk(maxloc(cnt),cnt);
        end
    end
    for cnt = 1:length(loc)
        pos = round((loc(cnt)-fqscale(1))/res)+1;
        T_f(pos,cnt) = val(cnt);
    end
    pic = figure;imagesc([0 timeEnd],[lowBound/res upBound/res],(T_f));colorbar off;title('SST');axis square;axis xy;
    xlabel('x'); ylabel('k');
    pic = figure;imagesc([0 timeEnd],[lowBound/res upBound/res],log2(T_f));colorbar off;title('SST');axis square;axis xy;
    xlabel('x'); ylabel('k');
end
