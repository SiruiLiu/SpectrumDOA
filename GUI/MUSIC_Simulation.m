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
    FFT_NUM = 256;
    FFT_NUM_LONG = 65536;
    d1 = figure(1001);
    d2 = figure(1002);
    d3 = figure(2000);
    d4 = figure(2001);
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
        for j = 1:element_num
           [S, F, T] = spectrogram(IQData(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, 9.6e6);
           figure(d1);
           subplot(2, element_num, j);
           imagesc(F, T, fftshift(db(abs(S.')), 2));
        end
        
        for j = 1:element_num
           [S, F, T] = spectrogram(CaliIQData(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, Fs);
           [S_IQ] = spectrogram(IQData(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, Fs);
           if(j==1)
               S1 = S;
               S_IQ1 = S_IQ;
           end
           figure(d1);
           subplot(2, element_num, j+element_num);
           imagesc(F, T, fftshift(db(abs(S.')), 2));
           if(j > 1)
                figure(d2)
                subplot(1, element_num-1, j-1);
                imagesc(F, T, fftshift(angle(S.'.*conj(S1.')),2));
                title('校正数据相位差');
                figure(d3)
                subplot(1, element_num-1, j-1);
                imagesc(F, T, fftshift(angle(S_IQ.'.*conj(S_IQ1.')),2));
                title('IQ数据相位差');
           end
        end
        
       %% Phase Calibration
        % Use the phase of maximum point as channel phase
        % For the current system, each ADC may exists an uncertain delay, this
        % time delay must be calibrated.
        % ---------------------------------------------------
        [CaliNew, IQDataNew] = DataCalibration(CaliIQData, IQData);
        
       %% Check if the calibration is correct
        for j = 1:element_num
           [S, F, T] = spectrogram(CaliNew(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, Fs);
           [S_IQ] = spectrogram(IQDataNew(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, Fs);
           if(j==1)
               S1 = S;
               S_IQ1 = S_IQ;
           end
           figure(d1);
           subplot(2, element_num, j+element_num);
           imagesc(F, T, fftshift(db(abs(S.')), 2));
           if(j > 1)
                figure(d2)
                subplot(1, element_num-1, j-1);
                imagesc(F, T, fftshift(angle(S.'.*conj(S1.')),2));
                title('延时校正后数据相位差');
                figure(d3)
                subplot(1, element_num-1, j-1);
                imagesc(F, T, fftshift(angle(S_IQ.'.*conj(S_IQ1.')),2));
                title('延时校正后IQ数据相位差');
           end
        end
                        
       %% Calculates phase difference.
        CaliFFTResult = fft(CaliNew, 32768, 2);
        CaliFreqSpec = db(abs(CaliFFTResult));
        [~, p] = max(CaliFreqSpec, [], 2);
        PhaseDiff = zeros(1, element_num - 1);
        for i = 1:element_num -1
            PhaseDiff(i) = angle(CaliFFTResult(i+1, p(i+1))*conj(CaliFFTResult(1, p(1))));
        end
           
       %% IQ Data filtering
        if(strcmp(Filter, 'AddFilter'))
            Coe1tmp = load('..\FilterCOE\Coe.mat'); % Loapass filter coefficients stored in .mat file
            coe1 = Coe1tmp.coe;
            IQDataFiltered = filter(coe1, 1, IQDataNew, [], 2);
            figure(1004);
            [S, F, T] = spectrogram(IQDataFiltered(1,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, 9.6e6);
            imagesc(F, T, fftshift(db(abs(S.')), 2));
        elseif(strcmp(Filter, 'CutOff'))
            IQFFTResult = fft(IQDataNew, FFT_NUM_LONG, 2);
            CutPoint1 = 0.5*length(IQFFTResult(1,:))+round((2.8e6+0.15e6)/Fs*length(IQFFTResult(1,:)));
            CutPoint2 = 0.5*length(IQFFTResult(1,:))+round((2.8e6-0.15e6)/Fs*length(IQFFTResult(1,:)));
            IQFFTResult(:,1:CutPoint2) = 0;
            IQFFTResult(:,CutPoint1:end) = 0;
            IQDataFiltered = ifft(IQFFTResult, FFT_NUM_LONG, 2);
            figure(1004);
            [S, F, T] = spectrogram(IQDataFiltered(1,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, 9.6e6);
            imagesc(F, T, fftshift(db(abs(S.')), 2));
        else
            IQDataFiltered = IQDataNew;
        end    
        
        IQData_Calibrated = PhaseCalibration(IQDataFiltered, PhaseDiff);
    else
        IQData_Calibrated = SampledData;
    end
    
    for j = 1:element_num
       [S_IQ] = spectrogram(IQData_Calibrated(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, Fs);
       if(j==1)
           S_IQ1 = S_IQ;
       end
       if(j > 1)
            figure(d4)
            subplot(1, element_num-1, j-1);
            imagesc(F, T, fftshift(angle(S_IQ.'.*conj(S_IQ1.')),2));
            title('IQ数据相位差');
       end
    end
    
    %% Source number estimation
    if(nargin <= 5)
        signal_num = SourceEst(IQData_Calibrated);
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
%     UU=O(:,1:element_num-(1));
%     Pmusic = zeros(length(theta),length(phi));
%     for i = 1:length(theta)
%         for j = 1:length(phi)
%             unit = [sin(theta(i))*cos(phi(j)),sin(theta(i))*sin(phi(j)),cos(theta(i))];
%             SteeringVector = exp(1i*2*pi/Lambda*Position*unit');
%             WW = SteeringVector'*(UU*UU')*SteeringVector;
%             Pmusic(i,j) = abs(1/WW);
%         end
%     end

    % Pmusic=10*log10(Pmusic/max(max(Pmusic)));

    figure(102);
    meshc(Pmusic);
    grid on
    
    figure(103)   % 三维视图
    meshc(sin(theta)'*cos(phi),sin(theta)'*sin(phi),Pmusic);
    
    [Theroy, Real] = PhaseDiff_Verify([90, 0], IQData_Calibrated, Lambda, Position);
end

function [Covariance] = Toeplitz(SigArray)
    ElementsNum = length(SigArray(:, 1));
    M_Array = zeros(1, ElementsNum);
    for i = 1:ElementsNum
         M_Array(i) = mean(SigArray(1,:).*conj(SigArray(i,:)));
    end
    R_Array = zeros(ElementsNum, ElementsNum);
    for i = 1:ElementsNum
        for j = 1:ElementsNum
            for k = 0:ElementsNum - 1
                if((j - i) == k)
                    R_Array(i, j) = M_Array(k+1);
                elseif((i - j) == k)
                    R_Array(i, j) = conj(M_Array(k+1));
                end
            end
        end
    end
    Covariance = R_Array;
end

function [Covariance] =  MMusic(OriginalCovariance)
    J = flip(eye(size(OriginalCovariance)));
    Covariance = OriginalCovariance + J*((OriginalCovariance)^-1*det(OriginalCovariance))*J;
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

function [Calibrated] = PhaseCalibration(OriginData, PhaseDiff)
    CaliArray = ones(size(OriginData,1)-1, size(OriginData, 2));
    CaliArray = CaliArray.*exp(1j*PhaseDiff).';
    Calibrated =  [OriginData(1,:); OriginData(2:end,:).*conj(CaliArray)];
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
function [CaliData_Calibrated, IQData_Calibrated] = DataCalibration(CaliData, IQData)
    global FFT_NUM_LONG;
    FFTResult = fft(CaliData, FFT_NUM_LONG, 2);
    Prod = FFTResult(2:end,:).*conj(FFTResult(1,:));
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
