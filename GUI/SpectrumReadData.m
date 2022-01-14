function [Frame] = SpectrumReadData(Path)
    TotalCh = 6; %共6通道
    FrameCh = 2; %每帧含2通道数据
    FFT_NUM = 1024;
    ShortDataPerFrame = 512*4; %按照short类型采集，512个DWORD相当于512*4个DBYTE
    
    fp = fopen(Path);

    DataTmp = fread(fp, [1,4], 'int32');
    DataFrame = struct([]);
    DataFrame(1).FrameHead = DataTmp(2); %帧头
    DataFrame(1).FrameCnt = DataTmp(1);  %帧计数 
    DataFrame(1).FrameLen = DataTmp(3);  %帧长，固定4072
    DataFrame(1).FrameFlg = DataTmp(4);  %帧标识
    
    fseek(fp,0,-1);
    data1 = fread(fp, 'uint8');
    DataTmp = reshape(data1, 8, [])';
    DataTmp = fliplr(DataTmp);
    
    % 时间戳
    FrameTime = struct([]);
    FrameTime(1).Nonosec = DataTmp(3,8) + ...
                           bitshift(DataTmp(3, 7), 8, 'uint32') + ...
                           bitshift(DataTmp(3, 6), 16, 'uint32') + ...
                           bitshift(DataTmp(3, 5), 24, 'uint32');%纳秒
    FrameTime(1).Sec = DataTmp(3, 4); %秒
    FrameTime(1).Minute = DataTmp(3, 3); %分
    FrameTime(1).Hour = DataTmp(3, 2); %时
    DataFrame(1).FrameTime = FrameTime;
    
    % 中心频率
    CentralFreq = DataTmp(4,8) + ...
                  bitshift(DataTmp(4, 7), 8, 'uint32') + ...
                  bitshift(DataTmp(4, 6), 16, 'uint32') + ...
                  bitshift(DataTmp(4, 5), 24, 'uint32');
    DataFrame(1).CentralFreq = CentralFreq;
    
    % IQ数据
    fseek(fp, 0, -1);
    data1 = fread(fp, 'short');% IQ数据每个采样点为16bit，即short类型，故按照short类型读取
    TotalFrame = length(data1)/ShortDataPerFrame;
    DataCh12 = [];
    DataCh34 = [];
    DataCh56 = [];
    %共1200帧，6通道，每帧2通道，故每两通道为1200/(6/2)=400帧
    for i = 1:TotalFrame/(TotalCh/FrameCh)
        %通道1~2数据起始位置为第6DWORD,一个DWORD为4个short类型，故通道1起始位置为每帧第4*5+1=21个short数据。
        %根据FPGA方面反馈，IQ数据共占用500个DWORD，故通道1数据相对帧起始位置偏移500*4+21-1=2020个点，从协
        %议上也能看出。
        %每两个通道按3帧循环，故每两个通道两组数据之间共间隔512DWOR*4*3=2048short*3=6144short个点的数据。
        DataCh12 = [DataCh12;data1((i-1)*6144+21:(i-1)*6144+2020)];
        %通道2~3相对于通道1~2偏移512WORD*4=2048个short类型数据。
        DataCh34 = [DataCh34;data1((i-1)*6144+2048+21:(i-1)*6144+2048+2020)];
        DataCh56 = [DataCh56;data1((i-1)*6144+4096+21:(i-1)*6144+4096+2020)];
    end 
    
    DataCh12_IQ = reshape(DataCh12, 2, []);
    DataCh34_IQ = reshape(DataCh34, 2, []);
    DataCh56_IQ = reshape(DataCh56, 2, []);
    DataCh12_IQ = DataCh12_IQ(1,:) + 1j*DataCh12_IQ(2,:);
    DataCh34_IQ = DataCh34_IQ(1,:) + 1j*DataCh34_IQ(2,:);
    DataCh56_IQ = DataCh56_IQ(1,:) + 1j*DataCh56_IQ(2,:);
    DataCh12_IQ = reshape(DataCh12_IQ, 2, []);
    DataCh34_IQ = reshape(DataCh34_IQ, 2, []);
    DataCh56_IQ = reshape(DataCh56_IQ, 2, []);
    Ch1_IQ = DataCh12_IQ(1,:);
    Ch2_IQ = DataCh12_IQ(2,:);
    Ch3_IQ = DataCh34_IQ(1,:);
    Ch4_IQ = DataCh34_IQ(2,:);
    Ch5_IQ = DataCh56_IQ(1,:);
    Ch6_IQ = DataCh56_IQ(2,:);
    
    IQData=zeros(TotalCh, length(Ch1_IQ));
%     IQData(1,:) = Ch1_IQ;
%     IQData(2,:) = Ch2_IQ;
%     IQData(3,:) = Ch3_IQ;
%     IQData(4,:) = Ch4_IQ;
%     IQData(5,:) = Ch5_IQ;
%     IQData(6,:) = Ch6_IQ;
    IQData(1,:) = Ch1_IQ - mean(Ch1_IQ);
    IQData(2,:) = Ch2_IQ - mean(Ch2_IQ);
    IQData(3,:) = Ch3_IQ - mean(Ch3_IQ);
    IQData(4,:) = Ch4_IQ - mean(Ch4_IQ);
    IQData(5,:) = Ch5_IQ - mean(Ch5_IQ);
    IQData(6,:) = Ch6_IQ - mean(Ch6_IQ);
    DataFrame(1).IQData=IQData;
    Frame = DataFrame;
    fclose(fp);
    
%     figure;
%     for i = 1:TotalCh
%        [S, F, T] = spectrogram(IQData(i,:), hamming(FFT_NUM), FFT_NUM*0.5, FFT_NUM, 153.6e6);
%        subplot(1, TotalCh, i);
%        imagesc(F, T, db(abs(S.')));
%     end

end

