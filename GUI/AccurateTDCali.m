clear;
close all;
clc;

%% ==============================================================
%
%                       Prepare data                                              
%
% ========================================================================
global FFT_SHORT_NUM;
global FFT_NUM; 
global FFT_LONG_NUM;
global FFT_LONGLONG_NUM;
FFT_SHORT_NUM =256;
FFT_NUM = 512;
FFT_LONG_NUM = 4096;
FFT_LONGLONG_NUM = 4096;
Fs = 9.6e6;
CaliLen_start = 1;
CaliLen_stop = 0.07;
IQData_Start = 0.5;
FreqDiff_Source = 3e6;
FreqInterval = 4e6;

% --Read data
load('CalibrationData.mat');
[CaliData, IQData] = Divide(DataStream, CaliLen_start, ...
                            round(size(DataStream, 2)*CaliLen_stop), ...
                            round(size(DataStream, 2)*IQData_Start));
[CaliDataNew, IQDataNew] = DataCalibration(CaliData, IQData);            
% Paintting(DataStream, CaliDataNew, IQDataNew, Fs);
AlignedData = AccurateTDCalibration(CaliDataNew, Fs, 1000);
Paintting(DataStream, AlignedData, IQDataNew, Fs);
% AlignedData2 = AccurateTDCalibration2(CaliDataNew, Fs);
%% ==============================================================
%
%                       Functions declared below                                              
%
% ========================================================================
% --Devide data stream into calibration data and iq data
function [CaliData, IQData] = Divide(DataStream, CaliStart, CaliStop, IQStart)
    CaliData = DataStream(:, CaliStart:CaliStop);
    IQData = DataStream(:, IQStart:end);
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

function Paintting2(Data, FFTN, Fs)
    for j = 1:size(Data, 1)
        [S, F, T] = spectrogram(Data(j,:), hamming(FFTN), FFTN*0.5, FFTN, Fs);
        if(j==1)
           S1 = S;
        end
        if(j > 1)
            figure(50000)
            subplot(1, size(Data, 1)-1, j-1);
            imagesc(F, T, fftshift(angle(S.'.*conj(S1.')),2));
            title(['通道', num2str(j), '与通道1延时校正后校正数据相位差']);
        end
    end
end

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

function [AlignedData] = AccurateTDCalibration(CaliData, Fs, Mul)
    global FFT_NUM;
    FFTResult = fft(CaliData, FFT_NUM, 2);
    FFTHighRes = fft(CaliData, [], 2);
    PhaseDiffperChannel = zeros(size(CaliData, 1)-1, FFT_NUM);
    PhaseDiffArray = zeros(size(CaliData, 1)-1, FFT_NUM);
    figure(10000);
    weights = zeros(size(CaliData, 1)-1,FFT_NUM);
    for j = 1:size(CaliData, 1)
       [S, F, T] = spectrogram(CaliData(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, Fs);
       if(j==1)
           S1 = S;
       end
       if(j > 1)
%             PhaseDiff = fftshift(angle(S.'.*conj(S1.')),2);
            subplot(3, size(CaliData, 1)-1, j-1);
            PhaseDiffTmp = S.'.*conj(S1.');
            PhaseDiffForPaintting = angle(PhaseDiffTmp);
            imagesc(F, T, PhaseDiffForPaintting);
            title(['通道', num2str(j), '与通道1延时校正后校正数据相位差']);
            PhaseDiffArray(j-1,:) = sum(PhaseDiffTmp);
            PhaseDiffperChannel(j-1,:) = angle(PhaseDiffArray(j-1,:));
            subplot(3, size(CaliData, 1)-1, j+size(CaliData, 1)-2);
            plot(PhaseDiffperChannel(j-1,:));
            subplot(3, size(CaliData, 1)-1, j+2*(size(CaliData, 1)-1)-1);
            plot(db(abs(FFTResult(j,:))));
            std_along_time = std(exp(1j*PhaseDiffForPaintting),[],1);
            std_along_time = -(std_along_time-max(std_along_time))/(max(std_along_time)-min(std_along_time))*5;
            weights(j-1,:) = exp(std_along_time)-exp(min(std_along_time));
       end
    end   
 
    PhaseDiffArray = PhaseDiffArray./abs(PhaseDiffArray);
    PhaseDiffArray(:,2:275) = 0;
    gcc = ifft(PhaseDiffArray, Mul*size(PhaseDiffArray, 2), 2);
    gcc_amp = abs(gcc);
    [~, p] = max(gcc_amp, [], 2);
    delay_point = zeros(length(p), 1);
    N = size(PhaseDiffArray, 2);
    omega = fftshift(2*pi/N*(-(N/2):(N/2-1)));
    for i = 1:length(p)
        if((p(i)>0) && (p(i)<size(gcc_amp, 2)*0.5))
            delay_point(i)=p(i)-1;
        else
            delay_point(i)=p(i)-1-size(gcc_amp, 2);
        end
        subplot(3,7,7+i),hold on
        calib_phase = -omega*delay_point(i)/Mul+angle(PhaseDiffArray(i,1));
        plot(calib_phase,'r');
    end
    Tau = 1/(Mul*Fs).*delay_point;
    AlignedData = Calibration(FFTHighRes, -1.*Tau, Fs);
end

function [AlignedData] = Calibration(CaliData_FFTResult, Tau, Fs)
    F = fftshift(linspace(Fs*0.5-Fs, Fs*0.5-Fs/size(CaliData_FFTResult, 2), size(CaliData_FFTResult, 2)));
    Omega = exp(1j*2*pi.*F.*Tau);
    FFT_Calibrated = zeros(size(CaliData_FFTResult, 1),size(CaliData_FFTResult, 2));
    FFT_Calibrated(1,:) = CaliData_FFTResult(1,:);
    FFT_Calibrated(2:end,:) = CaliData_FFTResult(2:end, :).*conj(Omega);
    AlignedData = ifft(FFT_Calibrated, [], 2);
end