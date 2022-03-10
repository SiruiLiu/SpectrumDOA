function [signal_num] = MUSIC_Simulation(varargin)
    %% 获取函数参数
    SampledData = varargin{1};
    Position = varargin{2};
    Lambda = varargin{3};
    Flg = 0;
    Filter = 0;
           
    if(nargin == 4)
        Flg = varargin{4};
    elseif(nargin == 5)
        Filter = varargin{5};
        Flg = varargin{4};
    elseif(nargin == 6)
        signal_num = varargin{6};
        Flg = varargin{4};
        Filter = varargin{5};
    end
    
    %% 参数设置
    global FFT_NUM;
    global FFT_NUM_LONG;
    FFT_NUM = 256;
    FFT_NUM_LONG = 65536;
    FreqDiff = 4.7e6;
    Fs = 9.6e6;
    element_num = length(SampledData(:,1));
    %% 切分数据 
    CaliLen_start = size(SampledData, 2)*0+1;
    CaliLen_stop = size(SampledData, 2)*0.07;
    IQData_Start = size(SampledData, 2)*0.5;

    angle_step = pi/180;       % 角度搜索步进 
    theta = angle_step:angle_step:pi/2;     % 俯仰角搜索范围
    phi = angle_step:angle_step:2*pi;     % 方位角搜索范围       
    
    %% Use real sampled data
    if(strcmp(Flg, 'REAL'))
        % Acquire calibration data and IQ data
        CaliIQData = SampledData(:, CaliLen_start:CaliLen_stop);
        IQData = SampledData(:, IQData_Start:end);
       %% Phase Calibration
        % Use the phase of maximum point as channel phase
        % For the current system, each ADC may exists an uncertain delay, this
        % time delay must be calibrated.
        % ---------------------------------------------------
        [CaliDataNew, IQDataNew] = DataCalibration(CaliIQData, IQData);
        [CaliData_Aligned, IQData_Aligned] = AccurateTDCalibration(CaliDataNew, IQDataNew, Fs, 100);
        Paintting(SampledData, CaliData_Aligned, IQData_Aligned, Fs);
        IQData_Calibrated = CaliData_Aligned;           
    else
        IQData_Calibrated = SampledData;
    end
    
    %% Source number estimation
    if(nargin <= 5)
        signal_num = SourceEst(IQData_Aligned);
    end
    
    %% MUSIC Core algorithm
    Covariance = IQData_Calibrated*IQData_Calibrated'./size(IQData_Calibrated, 2);
    [O, ~] = eig(Covariance); 
    O = O(1:element_num, 1:element_num);        

    % 谱峰搜索
    UU=O(:,1:element_num-(signal_num));
    Pmusic = zeros(length(theta),length(phi));
    for i = 1:length(theta)
        for j = 1:length(phi)
            unit = [sin(theta(i))*cos(phi(j)),sin(theta(i))*sin(phi(j)),cos(theta(i))];
            SteeringVector = exp(-1i*2*pi/Lambda*Position*unit');
            WW = SteeringVector'*(UU*UU')*SteeringVector;
            Pmusic(i,j) = abs(1/WW);
        end
    end

    figure(102);
    meshc(Pmusic);
    grid on
    
    figure(103)   % 三维视图
    meshc(sin(theta)'*cos(phi),sin(theta)'*sin(phi),Pmusic);
    
%     [Theroy, Real] = PhaseDiff_Verify([90, 0], IQData_Calibrated, Lambda, Position);
end

function [PhaseDiff] = CalPhaseDiff(IQData)
    global FFT_NUM_LONG;
    FFTResult = fft(IQData, FFT_NUM_LONG, 2);
    ampFFTResult = abs(FFTResult);
    [~, p] = max(ampFFTResult, [], 2);
    PhaseDiff = zeros(1, size(IQData, 1));
    for i = 1:size(IQData, 1)-1
        PhaseDiff(i) = angle(FFTResult(i+1, p(i+1))*conj(FFTResult(1, p(1))));
    end
end

function [PhaseDiff_Theory, PhaseDiff_RealData] = PhaseDiff_Verify(Ang, RealData, Lambda, Position)
    Radius = Ang*pi/180;
    unit = [sin(Radius(1))*cos(Radius(2)),sin(Radius(1))*sin(Radius(2)),cos(Radius(1))];
    SteeringVector = exp(1i*2*pi/Lambda*Position*unit');
    PhaseDiff_Theory = zeros(1, size(RealData, 1) - 1);
    for i = 1:size(RealData, 1)-1
        PhaseDiff_Theory(i) = angle(SteeringVector(i+1)*conj(SteeringVector(1)));
    end
    PhaseDiff_RealData = CalPhaseDiff(RealData);
end

% --Devide data stream into calibration data and iq data
function [CaliData, IQData] = Devide(DataStream, CaliStart, CaliStop, IQStart)
    CaliData = DataStream(:, CaliStart:CaliStop);
    IQData = DataStream(:, IQStart:end);
end

% --Calibration time delay for both calibration data and iq data
% --Calibration time delay for both calibration data and iq data
function [CaliData_Calibrated, IQData_Calibrated] = DataCalibration(CaliData, IQData)
    global FFT_NUM_LONG;
    FFTResult = fft(CaliData, FFT_NUM_LONG, 2);
    amp_tmp = abs(FFTResult);
    phase = angle(FFTResult);
    amp = log(amp_tmp)-min(log(amp_tmp), [], 2);
    ResultNew = amp.*exp(1j*phase);
    Prod = ResultNew(2:end,:).*conj(ResultNew(1,:));
    gcc = ifft(Prod, FFT_NUM_LONG, 2);
    gcc_amp = abs(gcc);

    [~, p] = max(gcc_amp, [], 2);
    delay_point = zeros(length(p), 1);
    for i = 1:length(p)
        if((p(i)>0) && (p(i)<FFT_NUM_LONG*0.5))
            delay_point(i)=p(i)-1;
        else
            delay_point(i)=p(i)-1-FFT_NUM_LONG;
        end
    end
    CaliData_Calibrated = TimeDelayCalibration(CaliData, delay_point);
    IQData_Calibrated = TimeDelayCalibration(IQData, delay_point);
end

% --Time delay calibration function
function [DataCalibrated] = TimeDelayCalibration(Data, delay_point)
    if(size(Data, 1)-1 ~= length(delay_point))
        error('Channel number and delay point number are mismatch');
    end
    CutPoint = max(delay_point) - min(delay_point);
    LagPoint = abs(min(delay_point));
    AheadPoint = max(delay_point);
    DataCalibrated = zeros(size(Data, 1), length(Data(1,:))-CutPoint);
    DataCalibrated(1,:) = Data(1, 1+LagPoint:end-AheadPoint);
    for i = 1:length(delay_point)
        if(delay_point(i) == 0)
            DataCalibrated(i+1,:) = Data(i+1, 1+LagPoint:end-AheadPoint);
        elseif(delay_point(i) < 0)
            DataCalibrated(i+1,:) = Data(i+1, 1:end-CutPoint);
        elseif(delay_point(i) > 0)
            DataCalibrated(i+1,:) = Data(i+1, 1+CutPoint:end);
        end
    end
end

% --Paintting function
function Paintting(DataStream, CaliData, IQData, Fs)
    global FFT_NUM;
    element_num = size(DataStream, 1);
    Diag1 = figure(1);
    Diag2 = figure(2);
    Diag3 = figure(3);
    
    for j = 1:element_num
       [S, F, T] = spectrogram(DataStream(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, 9.6e6);
       figure(Diag1);
       subplot(1, element_num, j);
       imagesc(F, T, fftshift(db(abs(S.')), 2));
       title(['通道', num2str(j),'数据时频图']);
    end
    
    for j = 1:element_num
       [S] = spectrogram(CaliData(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, Fs);
       [S_IQ] = spectrogram(IQData(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, Fs);
       figure(Diag2);
       subplot(2, element_num, j);
       imagesc(F, T, fftshift(db(abs(S.')), 2));
       title(['通道', num2str(j), '校正数据时频图']);
       subplot(2, element_num, j+element_num);
       imagesc(F, T, fftshift(db(abs(S_IQ.')), 2));
       hold on;
       [v] = max(max(fftshift(db(abs(S_IQ.')))));
       [x, y] = find(fftshift(db(abs(S_IQ.'))) == v);
       plot(F(y), T(x), 'r+');
       text(F(y), T(x), num2str(v))
       hold off;
       title(['通道', num2str(j), 'IQ数据时频图']);
       if(j==1)
           S1 = S;
           S_IQ1 = S_IQ;
       end
       if(j > 1)
            figure(Diag3)
            subplot(2, element_num-1, j-1);
            imagesc(F, T, fftshift(angle(S.'.*conj(S1.')),2));
            title(['通道', num2str(j), '与通道1延时校正后校正数据相位差']);
            subplot(2, element_num-1, j+element_num-2);
            imagesc(F, T, fftshift(angle(S_IQ.'.*conj(S_IQ1.')),2));
            title(['通道', num2str(j), '与通道1延时校正后IQ数据相位差']);
       end
    end   
end

function [CaliData_AlignedData, IQData_AlignedData] = AccurateTDCalibration(CaliData, IQData, Fs, Mul)
    global FFT_NUM;
    CaliData_FFTHighRes = fft(CaliData, [], 2);
    IQData_FFTHighRes = fft(IQData, [], 2);
    PhaseDiffperChannel = zeros(size(CaliData, 1)-1, FFT_NUM);
    PhaseDiffArray = zeros(size(CaliData, 1)-1, FFT_NUM);
    for j = 1:size(CaliData, 1)
       [S] = spectrogram(CaliData(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, Fs);
       if(j==1)
           S1 = S;
       end
       if(j > 1)
            PhaseDiffTmp = S.'.*conj(S1.');
            PhaseDiffArray(j-1,:) = sum(PhaseDiffTmp);
            PhaseDiffperChannel(j-1,:) = angle(PhaseDiffArray(j-1,:));
       end
    end   
 
    PhaseDiffArray = PhaseDiffArray./abs(PhaseDiffArray);
    PhaseDiffArray(:,2:275) = 0;
    gcc = ifft(PhaseDiffArray, Mul*size(PhaseDiffArray, 2), 2);
    gcc_amp = abs(gcc);
    [~, p] = max(gcc_amp, [], 2);
    delay_point = zeros(length(p), 1);
    N = size(PhaseDiffArray, 2);
    for i = 1:length(p)
        if((p(i)>0) && (p(i)<size(gcc_amp, 2)*0.5))
            delay_point(i)=p(i)-1;
        else
            delay_point(i)=p(i)-1-size(gcc_amp, 2);
        end
    end
    Tau = 1/(Mul*Fs).*delay_point;
    CaliData_AlignedData = Calibration(CaliData_FFTHighRes, -1.*Tau, Fs);
    IQData_AlignedData = Calibration(IQData_FFTHighRes, -1.*Tau, Fs);
end

function [AlignedData] = Calibration(CaliData_FFTResult, Tau, Fs)
    F = fftshift(linspace(Fs*0.5-Fs, Fs*0.5-Fs/size(CaliData_FFTResult, 2), size(CaliData_FFTResult, 2)));
    Omega = exp(1j*2*pi.*F.*Tau);
    FFT_Calibrated = zeros(size(CaliData_FFTResult, 1),size(CaliData_FFTResult, 2));
    FFT_Calibrated(1,:) = CaliData_FFTResult(1,:);
    FFT_Calibrated(2:end,:) = CaliData_FFTResult(2:end, :).*conj(Omega);
    AlignedData = ifft(FFT_Calibrated, [], 2);
end