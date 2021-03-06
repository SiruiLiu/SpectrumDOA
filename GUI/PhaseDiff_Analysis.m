clear;
close all;
clc;

%% =========Parameters declared below=======
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

PositionStd=[0.15,        0, 0;
             0.106,  -0.106, 0;
             0,       -0.15, 0;
             -0.106, -0.106, 0;
             -0.15,       0, 0;
             -0.106,  0.106, 0;
             0,        0.15, 0;
             0.106,   0.106, 0];
loc_theta = -45*(0:1:7);
loc_theta = loc_theta(:);
r = 0.15;
PositionStd=[r*cosd(loc_theta),r*sind(loc_theta),zeros(8,1)];

ChannelNum = 8;
Position = PositionStd(1:ChannelNum, :);
Height = 1.1;
Dist = 16;
% Angle = [atan(Dist/Height)*180/pi, 360-135];
Ang  = [90, 0;
        90, -15;
        90, -30;
        90, -45;
        90, -60;
        90, -75;
        90, -90;];

% Calibration data and IQData a construct in one data stream, use the 
% parameters below to devide them
CaliLen_start = 1;
CaliLen_stop = 0.07;
IQData_Start = 0.5;

Diag4 = figure(1);
Diag5 = figure(2);
Diag6 = figure(3);
%% ================Visit data file directory===========
FilePath = uigetdir(' F:\SpectrumDOA\SimData\', 'Select a directory');
FileInfo = dir(fullfile(FilePath));
FileInfo = FileInfo(3:end);
% FileInfo = FileInfo(1);
% FileInfo = FileInfo(1:end-2);
FileList = {FileInfo.name};
%% ===============Comparing theory phase difference and real phase difference============
Diff = zeros(length(FileList), ChannelNum-1);
k = 1;
AngleStream = [];
for i = 1:length(FileInfo)
    Path_2nd = [FileInfo(i).folder, '\', FileInfo(i).name];
    FileInfo_2nd = dir(fullfile(Path_2nd, '*.dat'));
    Angle = Ang(1,:);
    AngleStream = [AngleStream, abs(Ang(i, 2))];
    for j = 1:length(FileInfo_2nd)
%     for j = 22:23
        Path = [FileInfo_2nd(j).folder, '\', FileInfo_2nd(j).name];
        Frame = SpectrumReadData(Path);
        [CaliData, IQData] = Devide(Frame.IQData, CaliLen_start, ...
                                    round(size(Frame.IQData, 2)*CaliLen_stop), ...
                                    round(size(Frame.IQData, 2)*IQData_Start));
%         [CaliDataNew, IQDataNew] = DataCalibration(CaliData, IQData);
%         [CaliData_Aligned, IQData_Aligned] = AccurateTDCalibration(CaliDataNew, IQDataNew, Fs, 100);
        Paintting(Frame.IQData, CaliData, IQData, Fs);
        IQPhase = CalPhase(IQData);
%         CaliPhaseDiff = Cal_CaliPhaseDiff(CaliData_Aligned);
%         Paintting(Frame.IQData, CaliData_Aligned, IQData_Aligned, Fs);
%         IQPhase = CalPhase(IQDataNew);
%         CaliPhaseDiff = Cal_CaliPhaseDiff(CaliDataNew);
%         SteeringVec = CalTheoryPhase(Angle, Position);
%         Tmp = IQPhase.*conj(CaliPhaseDiff);
%         PhaseDiff(:,k) = angle(Tmp.*conj(Tmp(3)));
%         CaliPhaseDiff_Radius(:,k) = angle(CaliPhaseDiff);
%         k = k+1;
%         PhaseDiff(:,k) = unwrap(CalResidual(IQDataNew, CaliDataNew, Angle));
%         k = k+1;
    end
end
% PhaseDiff = unwrap(PhaseDiff, [], 2);
% CaliPhaseDiff_Radius = unwrap(CaliPhaseDiff_Radius, [], 2);
Variance = std(PhaseDiff*180/pi);
figure(Diag4);
for i = 1:size(PhaseDiff, 1)
    plot(PhaseDiff(i,:));
    hold on;
end
set(gca, 'xticklabel', [0, 15, 30, 45, 60, 75 ,90]);
legend('??????1?????????1?????????','??????2?????????1?????????','??????3?????????1?????????', ...
       '??????4?????????1?????????','??????5?????????1?????????','??????6?????????1?????????');
hold off;

figure(Diag5);
for i = 1:size(CaliPhaseDiff_Radius, 1)
    plot(CaliPhaseDiff_Radius(i,:));
    hold on;
end
set(gca, 'xticklabel', [0, 15, 30, 45, 60, 75 ,90]);
legend('??????1?????????1???????????????','??????2?????????1???????????????','??????3?????????1???????????????', ...
       '??????4?????????1???????????????','??????5?????????1???????????????','??????6?????????1???????????????');
hold off;

% figure(Diag5);
% for i = 1:size(PhaseDiff, 1)
%     plot(PhaseTmp1(i,:));
%     hold on;
% end
% set(gca, 'xticklabel', [0, 15, 30, 45, 60, 75 ,90]);
% legend('FFT??????','??????2','??????3', ...
%        '??????4','??????5','??????6');
% title('FFT??????');
% hold off;;
%% ========================Functions declared and described below =====================
% --Calculates theroy and real phase difference between channel 1 and
% --other channels
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

% --Calculates real data phase difference
function [PhaseDiff] = CalPhaseDiff(IQData)
    global FFT_NUM_LONG;
    FFTResult = fft(IQData, FFT_NUM_LONG, 2);
    ampFFTResult = abs(FFTResult);
    [~, p] = max(ampFFTResult, [], 2);
    PhaseDiff = zeros(1, size(IQData, 1)-1);
    for i = 1:size(IQData, 1)-1
        PhaseDiff(i) = angle(FFTResult(i+1, p(i+1))*conj(FFTResult(1, p(1))));
    end
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

% -- Phase calibration function
function [Calibrated] = PhaseCalibration(OriginData, PhaseDiff)
    CaliArray = ones(size(OriginData,1)-1, size(OriginData, 2));
    CaliArray = CaliArray.*exp(1j*PhaseDiff).';
    Calibrated =  [OriginData(1,:); OriginData(2:end,:).*conj(CaliArray)];
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
       title(['??????', num2str(j),'???????????????']);
    end
    
    for j = 1:element_num
       [S] = spectrogram(CaliData(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, Fs);
       [S_IQ] = spectrogram(IQData(j,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, Fs);
       figure(Diag2);
       subplot(2, element_num, j);
       imagesc(F, T, fftshift(db(abs(S.')), 2));
       title(['??????', num2str(j), '?????????????????????']);
       subplot(2, element_num, j+element_num);
       imagesc(F, T, fftshift(db(abs(S_IQ.')), 2));
       title(['??????', num2str(j), 'IQ???????????????']);
       if(j==1)
           S1 = S;
           S_IQ1 = S_IQ;
       end
       if(j > 1)
            figure(Diag3)
            subplot(2, element_num-1, j-1);
            imagesc(F, T, fftshift(angle(S.'.*conj(S1.')),2));
            title(['??????', num2str(j), '?????????1????????????????????????????????????']);
            subplot(2, element_num-1, j+element_num-2);
            imagesc(F, T, fftshift(angle(S_IQ.'.*conj(S_IQ1.')),2));
            title(['??????', num2str(j), '?????????1???????????????IQ???????????????']);
       end
    end   
end

% --Calculates steering vector(or the theory phase)
function [TheoryPhase] = CalTheoryPhase(Ang, Position)
    global Lambda;
    Radius = Ang*pi/180;
    unit = [sin(Radius(1))*cos(Radius(2)),sin(Radius(1))*sin(Radius(2)),cos(Radius(1))];
    TheoryPhase = exp(1i*2*pi/Lambda*Position*unit');
end

% --Use one radiation source signal to calibrating the other one, and check
% --if the phase differences are equal for difference direction.
% function [PhaseDiff1, PhaseDiff2] = FindPhaseDiff_for_2Freq(varargin)
%     global FFT_NUM_LONG;
%     global Position;
%     Data = varargin{1};
%     Fs = varargin{2};
%     Res = Fs/FFT_NUM_LONG;
%     Fc = varargin{3};
%     Freq1Info = varargin{4};
%     Freq2Info = varargin{5};
%     
%     DataFFTResult = fft(Data, FFT_NUM_LONG, 2);
%     ampFFTResult = db(abs(DataFFTResult));
%     [~, Point] = findpeaks(ampFFTResult(1,:), 'SortStr','descend','NPeaks', 2);
%     Point = sort(Point);
%     
%     Ang1 = Freq1Info(2:3);
%     Ang2 = Freq2Info(2:3);
%     
%     TheoryPhase1 = CalTheoryPhase(Ang1, Position);
%     TheoryPhase2 = CalTheoryPhase(Ang2, Position);
%     tmp1 = angle(TheoryPhase1);
%     tmp2 = angle(TheoryPhase2);
%     PhaseArray1 = DataFFTResult(:, Point(1));
%     PhaseArray2 = DataFFTResult(:, Point(2));
%     
%     PhaseDiff1_tmp = PhaseArray1.*conj(TheoryPhase1);
%     PhaseDiff2_tmp = PhaseArray2.*conj(TheoryPhase2);
%     PhaseDiff1 = angle(PhaseDiff1_tmp*conj(PhaseDiff1_tmp(1)));
%     PhaseDiff2 = angle(PhaseDiff2_tmp*conj(PhaseDiff2_tmp(1)));
% end

% --Phase calculation, assume the maximum power point is the signal's
% --frequency
function [Phase] = CalPhase(Data)
    global FFT_NUM_LONG;
    FFTResult = fft(Data, FFT_NUM_LONG, 2);
    ampFFTResult = abs(FFTResult);
    [~, p] = max(ampFFTResult, [], 2);
    Phase = zeros(size(Data, 1), 1);
    for i = 1:length(Phase)
        Phase(i) = FFTResult(i, p(i));
    end
end

% --
% function [PhaseDiff] = CalResidual(IQData, CaliData, Ang)
%     global Position;
%     SteeringVec = CalTheoryPhase(Ang, Position);
%     IQPhase = CalPhase(IQData);
%     CaliPhase = Cal_CaliPhase1(CaliData);
%     PhaseDiff = angle(IQPhase.*conj(SteeringVec).*conj(CaliPhase));
%     PhaseDiff = mod(PhaseDiff-PhaseDiff(1),2*pi);
% end

% function [CaliPhase] = Cal_CaliPhase1(CaliData)
%     global FFT_NUM_LONG;
%     FFTResult = fft(CaliData, FFT_NUM_LONG, 2);
%     PhaseDiff_Tmp = FFTResult.*conj(FFTResult(1,:));
%     CaliPhaseDiff = mean(PhaseDiff_Tmp.*abs(FFTResult), 2);
% end

function [CaliPhaseDiff] = Cal_CaliPhaseDiff(CaliData)
    global FFT_NUM_LONG;
    FFTResult = fft(CaliData, FFT_NUM_LONG, 2);
    PhaseDiff_Tmp = FFTResult.*conj(FFTResult(3,:));
    CaliPhaseDiff = mean(PhaseDiff_Tmp.*abs(FFTResult), 2);
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