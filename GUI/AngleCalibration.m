function [Result] = AngleCalibration(varargin)
%% -------Function for space spectrum calibration------
% This function used for mechanical error calibration, includes
% antenna's position calibration and input direction.
%% ------------Parameters preparing------------
    DirInfo = varargin{1};
    BaseCh = varargin{2};
    % Calibration data and IQData a construct in one data stream, use the 
    % parameters below to devide them
    CaliLen_start = varargin{3};
    CaliLen_stop = varargin{4};
    IQData_Start = varargin{5};
    Fs = varargin{6};
    Mul = varargin{7};
    Position = varargin{8};
    OutliersTh = varargin{9};
    
    global FFT_NUM;
    global C;
    FFT_NUM = 512;
    C=299792458;
%% ---------------------------------------------------------- %%
%
%                           Main Program
%
% -------------------------------------------------------------------- %%
    FileName = natsort({DirInfo.name});
    RelativeCoordinates = CalRelativeCoordinates(Position, BaseCh); %Calculates relative coordinate
    CoeArray_Z = zeros(length(DirInfo), size(Position, 1)-1);
    CoeArray_A1 = zeros(length(DirInfo), size(Position, 1)-1);
    CoeArray_A2 = zeros(length(DirInfo), size(Position, 1)-1);
    CoeArray_A3 = zeros(length(DirInfo), size(Position, 1)-1);
    CoeArray_B = zeros(length(DirInfo), size(Position, 1)-1);
    CoeArray_C = zeros(length(DirInfo), size(Position, 1)-1);
    PhaseDiff = [];
    IQPhaseDiff = [];
    AngleTmp = [];
    for i = 1:length(DirInfo)
        FilePath = [DirInfo(i).folder, '\', FileName{i}];
        FileInfo = dir(fullfile(FilePath, '*.dat'));
        Pos = strfind(FileName{i}, 'D');
        Angle = [0-str2double(FileName{i}(1:Pos-1)), 90];%俯仰角默认为90度。
        [TmpPhase, TmpIQPhase, AngleStream] = BuildFunction(FileInfo, Angle, BaseCh, ...
                                                CaliLen_start, CaliLen_stop, IQData_Start, ...
                                                Fs, Mul, RelativeCoordinates, OutliersTh);
        PhaseDiff = [PhaseDiff, TmpPhase];
        IQPhaseDiff = [IQPhaseDiff, TmpIQPhase];
        AngleTmp = [AngleTmp, AngleStream];
    end
    figure(1000);
    for i = 1:size(PhaseDiff, 1)
        subplot(1,4,i);
        plot(AngleTmp, PhaseDiff(i,:),'r.');
        title(['通道',num2str(i)]);
    end
    figure(1001);
    for i = 1:size(PhaseDiff, 1)
        subplot(1,4,i);
        plot(AngleTmp, IQPhaseDiff(i,:),'r.');
        title(['通道',num2str(i)]);
    end
    
%         [Z, CoeA1, CoeA2, CoeA3, CoeB, CoeC] = BuildFunction(FileInfo, Angle, BaseCh, ...
%                                                              CaliLen_start, CaliLen_stop, IQData_Start, ...
%                                                              Fs, Mul, RelativeCoordinates, OutliersTh);
                                                         
%         CoeArray_Z(i,:) = Z.';
%         CoeArray_A1(i,:) = repmat(CoeA1, size(Position, 1)-1, 1).';
%         CoeArray_A2(i,:) = repmat(CoeA2, size(Position, 1)-1, 1).';
%         CoeArray_A3(i,:) = repmat(CoeA3, size(Position, 1)-1, 1).';
%         CoeArray_B(i,:) = CoeB;
%         CoeArray_C(i,:) = CoeC;
%     end
%     Result = zeros(6, size(Position, 1)-1);
%     for i = 1:size(Position, 1)-1
%         CoeArray = [CoeArray_A1(:, i), CoeArray_A2(:, i), CoeArray_A3(:,i), CoeArray_B(:,i), CoeArray_C(:,i), ones(length(DirInfo),1)];
%         Result(:, i) = ((CoeArray.'*CoeArray)^-1)*CoeArray.'*CoeArray_Z(:,i);  
end

%% ---------------------------------------------------------- %%
%
%                      Functions declared below
%
% -------------------------------------------------------------------- %%
function [RelativeCoordinates] = CalRelativeCoordinates(Position, BaseCh)
    BasePosition = Position(BaseCh,:);
    RelativeCoordinates = Position;
    RelativeCoordinates(BaseCh,:)=[];
    BasePosition = repmat(BasePosition, size(RelativeCoordinates, 1), 1);
    RelativeCoordinates = RelativeCoordinates-BasePosition;
end

% function [Z, CoeA1, CoeA2, CoeA3, CoeB, CoeC] = BuildFunction(FileInfo, Angle, BaseCh, ...
%                                                 CaliLen_start, CaliLen_stop, IQData_Start, ...
%                                                 Fs, Mul, Position, OutliersTh)           
function [Phase_without_CaliPhase, IQPhase_without_CaliPhase, AngleStream] = BuildFunction(FileInfo, Angle, BaseCh, ...
                                                                                            CaliLen_start, CaliLen_stop, ...
                                                                                            IQData_Start, ...
                                                                                            Fs, Mul, Position, OutliersTh)     
    [Phase_without_CaliPhase, IQPhase_without_CaliPhase, AngleStream, Lambda] = PrepareData(FileInfo, BaseCh, ...
                                                                                            Angle, CaliLen_start, ...
                                                                                            CaliLen_stop, IQData_Start, ...
                                                                                            Fs, Mul, Position);
%     Z = CalStableMeanValue(Phase_without_CaliPhase, OutliersTh);
%     [CoeA1, CoeA2, CoeA3] = CalCoe_A(Angle, Lambda);
%     CoeB = CalCoe_B(Angle, Lambda, Position);
%     CoeC = CalCoe_C(Angle, Lambda, Position);
end


function [Phase_without_CaliPhase, IQPhase_without_CaliPhase, AngleStream, Lambda] = PrepareData(FileInfo, BaseCh, ...
                                                                      Angle, CaliStart, ...
                                                                      CaliStop, IQStart, ...
                                                                      Fs, Mul, Position)
    global C;
%     Phase_without_CaliPhase = zeros(size(Position, 1), length(FileInfo));
    Phase_without_CaliPhase = zeros(round(0.5*size(Position, 1)), length(FileInfo));
    IQPhase_without_CaliPhase = zeros(round(0.5*size(Position, 1)), length(FileInfo));
    AngleStream = zeros(1, length(FileInfo));
    for i = 1:length(FileInfo)
        Path = [FileInfo(i).folder, '\', FileInfo(i).name];
        Frame = SpectrumReadData(Path);
        Lambda = C/Frame.CentralFreq;
        [CaliData, IQData] =  Devide(Frame.IQData, CaliStart+1, ...
                                     round(size(Frame.IQData, 2)*CaliStop), ...
                                     round(size(Frame.IQData, 2)*IQStart));
        [CaliData_Aligned, IQData_Aligned] = CalibrationProcess(CaliData, IQData, ...
                                                                Fs, Mul);
%         Paintting(Frame.IQData, CaliData, IQData, Fs);
        IQPhase = CalPhase(IQData);
        AngleStream(i) = abs(Angle(1));
        for j = 1:round(0.5*size(Position, 1))
            IQPhase_without_CaliPhase(j, i) = angle(IQPhase(2*j).*conj(IQPhase(2*(j-1)+1)));
            Phase_without_CaliPhase(j, i) = abs(IQPhase(2*j).*conj(IQPhase(2*(j-1)+1)));
        end
%         CaliPhaseDiff = Cal_CaliPhaseDiff(CaliData_Aligned, BaseCh);
%         Phase_Calibrated = IQPhase.*conj(CaliPhaseDiff); %Remove ratio frequency channal phase difference
%         Phase_Calibrated = Phase_Calibrated.*conj(Phase_Calibrated(BaseCh));
%         Phase_Calibrated(BaseCh,:) = [];
%         SteeringVec = CalTheoryPhase(Angle, Position, Lambda);
%         Phase_without_CaliPhase(:, i) = angle(Phase_Calibrated.*conj(SteeringVec)); %Remove phase error casued by signal direction
%         Phase_without_CaliPhase(:, i) = angle(Phase_Calibrated);
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

function [CaliPhase] = CalCaliPhase(Data)
    global FFT_NUM_LONG;
    FFTResult = fft(Data, FFT_NUM_LONG, 2);
    CaliPhase = mean(FFTResult.*abs(FFTResult));
end

% --Calibration time delay for both calibration data and iq data
function [CaliData_Calibrated, IQData_Calibrated] = DataCalibration(CaliData, IQData)
    FFTResult = fft(CaliData, [], 2);
    amp_tmp = abs(FFTResult);
    phase = angle(FFTResult);
    amp = log(amp_tmp)-min(log(amp_tmp), [], 2);
    ResultNew = amp.*exp(1j*phase);
    Prod = ResultNew(2:end,:).*conj(ResultNew(1,:));
    gcc = ifft(Prod, size(CaliData, 2), 2);
    gcc_amp = abs(gcc);

    [~, p] = max(gcc_amp, [], 2);
    delay_point = zeros(length(p), 1);
    for i = 1:length(p)
        if((p(i)>0) && (p(i)<size(CaliData, 2)*0.5))
            delay_point(i)=p(i)-1;
        else
            delay_point(i)=p(i)-1-size(CaliData, 2);
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
function [CaliData_Aligned, IQData_Aligned] = CalibrationProcess(CaliData, IQData, ...
                                                                 Fs, Mul)
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

% --Calculates calibration data's phase difference
function [CaliPhaseDiff] = Cal_CaliPhaseDiff(CaliData, BaseCh)
    global FFT_NUM_LONG;
    FFTResult = fft(CaliData, FFT_NUM_LONG, 2);
    PhaseDiff_Tmp = FFTResult.*conj(FFTResult(BaseCh,:));
    CaliPhaseDiff = mean(PhaseDiff_Tmp.*abs(FFTResult), 2);
end

% --Calculates steering vector(or the theory phase)
function [TheoryPhase] = CalTheoryPhase(Ang, Position, Lambda)
    Radius = Ang*pi/180;
    unit = [sin(Radius(1))*cos(Radius(2)),sin(Radius(1))*sin(Radius(2)),cos(Radius(1))];
    TheoryPhase = exp(1i*2*pi/Lambda*Position*unit');
end

% --Due to the unstable phase, the point with large mistake should be
% --removed, and calculates the mean value.
function [Coe_Z] = CalStableMeanValue(Phase, Th)
    Standard = std(Phase, [], 2);
    R = find(Standard > Th);
    Coe_Z = mean(Phase, 2);
    Tmp = []; % Declared a empty matrix for adding legal value
    for i = 1:length(R)
        for j = 1:length(Phase(R(i),:))
            if((Phase(R(i),j) - mean(Phase(R(i),:))) < Standard(R(i)))
                Tmp = [Tmp,Phase(R(i),j)];
            end
        end
        Coe_Z(R(i))=mean(Tmp);
    end
end

function [A1, A2, A3] = CalCoe_A(Angle, Lambda)
    Angle = Angle*pi/180;
    A1 = 2*pi/Lambda*(sin(Angle(2))*cos(Angle(1)));
    A2 = 2*pi/Lambda*(sin(Angle(2))*sin(Angle(1)));
    A3 = 2*pi/Lambda*cos(Angle(2));
end

function [B] = CalCoe_B(Angle, Lambda, Position)
    Angle = Angle*pi/180;
    B = 2*pi/Lambda*([cos(Angle(2))*cos(Angle(1)), cos(Angle(2))*sin(Angle(1)), -sin(Angle(2))]*Position.');
end

function [C] = CalCoe_C(Angle, Lambda, Position)
    Angle = Angle*pi/180;
    C = 2*pi/Lambda*([-sin(Angle(2))*sin(Angle(1)), sin(Angle(2))*cos(Angle(1)), 0]*Position.');   
end
