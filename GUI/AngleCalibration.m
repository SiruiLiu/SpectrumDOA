function [Coe] = AngleCalibration(varargin)
%% -------Function for space spectrum calibration------
% This function used for mechanical error calibration, includes
% antenna's position calibration and input direction.
%% ------------Parameters preparing------------
    DirInfo = varargin{1};
    BaseCh = varargin{2};
    Angle = zeros(length(DirInfo), 2);
    % Calibration data and IQData a construct in one data stream, use the 
    % parameters below to devide them
    CaliLen_start = varargin{3};
    CaliLen_stop = varargin{4};
    IQData_Start = varargin{5};
    Fs = varargin{6};
    Mul = varargin{7};
%% ---------------------------------------------------------- %%
%
%                           Main Program
%
% -------------------------------------------------------------------- %%
    for i = 1:length(DirInfo)
        FilePath = [DirInfo(i).folder, '\', DirInfo(i).name];
        FileInfo = dir(fullfile(FilePath, '*.dat'));
        Pos = strfind(DirInfo(i).name, 'D');
        Angle = [num2str(FileInfo(i).name(1:Pos-1)), 90];%俯仰角默认为90度。
        IQData = BuildFunction(FileInfo, Angle(i,:), CaliLen_start, CaliLen_stop, IQData_Start, Fs, Mul);
    end
end

%% ---------------------------------------------------------- %%
%
%                      Functions declared below
%
% -------------------------------------------------------------------- %%
function [Coe_A, Coe_B] = BuildFunction(FileInfo, Angle, CaliLen_start, CaliLen_stop, IQData_Start, Fs, Mul)
    [CaliData, IQData] = PrepareData(FileInfo, CaliLen_start, CaliLen_stop, IQData_Start);
    [CaliData_Aligned, IQData_Aligned] = CalibrationProcess(CaliData, IQData, Fs, Mul);
    Paintting(IQData, CaliData_Aligned, IQData_Aligned, Fs);
end


function [CaliData, IQData] = PrepareData(FileInfo, CaliStart, CaliStop, IQStart)
    for i = 1:length(FileInfo)
        Path = [FileInfo(i).folder, '\', FileInfo(i).name];
        Frame = SpectrumReadData(Path);
        [CaliData, IQData] =  Devide(Frame.IQData(1,:), CaliStart+1, ...
                                                     round(size(Frame.IQData, 2)*CaliStop), ...
                                                     round(size(Frame.IQData, 2)*IQStart));
    end
end

% --Devide data stream into calibration data and iq data
function [CaliData, IQData] = Devide(DataStream, CaliStart, CaliStop, IQStart)
    CaliData = DataStream(:, CaliStart:CaliStop);
    IQData = DataStream(:, IQStart:end);
end

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

% --The preivous function could only calibrates the time delay longer than
% --sample period, this function could do further calibration for the delay
% --time short than sample period.
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

% --Do coarse calibration and accuracy calibration for both data
function [CaliData_Aligned, IQData_Aligned] = CalibrationProcess(CaliData, IQData, Fs, Mul)
    [CaliDataTmp, IQDataTmp] = DataCalibration(CaliData, IQData);
    [CaliData_Aligned, IQData_Aligned] = AccurateTDCalibration(CaliDataTmp, ...
                                                               IQDataTmp, ...
                                                               Fs, Mul);
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
