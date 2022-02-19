clear;
close all;
clc;

C = 299792458;
global FFT_NUM;
global FFT_NUM_LONG;
global Position;
global Lambda;
FFT_NUM = 512;
FFT_NUM_LONG = 65536;
Fc = 1103e6;
Fs = 9.6e6;
Freq1 = 1101e6;
Freq2 = 1102e6;
Lambda = C/Fc;

Diag4 = figure(1);
Diag5 = figure(2);
Diag6 = figure(3);

%% -----------Read Data----------
FilePath = uigetdir('..\sim_data\', 'Directory selector');
FileInfo = dir(fullfile(FilePath));
FileInfo = FileInfo(3:end);
% for i = 1:length(FileInfo)
%     p = strfind(FileInfo(i).name, '_');
%     FileIndex(i) = str2num(FileInfo(i).name(1:p(1)-1));
% end
% FileIndex = num2str(unique(FileIndex'));

%% ==============================================================
%
%                             MAIN PROCESS
%
%  =======================================================================
k = 1;
% for i = 1:size(FileIndex, 1)
%     Data = SDRDataRead(FilePath, FileIndex(i,:));
%     [CaliDataNew, IQDataNew] = DataCalibration(Data, Data);
%     Paintting(Data, Data, Data, Fs);
%     CaliPhaseDiff = Cal_CaliPhaseDiff(CaliDataNew);
%     CaliPhaseDiff_Radius(:,k) = unwrap(angle(CaliPhaseDiff));
%     k = k+1;
% end

for i = 1:length(FileInfo)
    File = [FileInfo(i).folder, '\', FileInfo(i).name];
    Data = SDRDataRead(File);
    [CaliDataNew, IQDataNew] = DataCalibration(Data, Data);
%     Paintting(Data, Data, Data, Fs);
    CaliPhaseDiff = Cal_CaliPhaseDiff(Data);
    CaliPhaseDiff_Radius(:,k) = unwrap(angle(CaliPhaseDiff));
    k = k+1;
end

figure(Diag6);
for i = 1:size(CaliPhaseDiff_Radius, 1)
    plot(CaliPhaseDiff_Radius(i,:));
    hold on;
end
legend('通道1与通道1校正相位差','通道2与通道1校正相位差','通道3与通道1校正相位差', ...
       '通道4与通道1校正相位差');
hold off;
%% ===============================================================
%
%                   FUNCTIONS DECLARED BELOW
%
% =========================================================================

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

function [CaliPhaseDiff] = Cal_CaliPhaseDiff(CaliData)
    global FFT_NUM_LONG;
    FFTResult = fft(CaliData, FFT_NUM_LONG, 2);
    PhaseDiff_Tmp = FFTResult.*conj(FFTResult(3,:));
    CaliPhaseDiff = mean(PhaseDiff_Tmp.*abs(FFTResult), 2);
end