% function [Data] = SDRDataRead(FilePath, FileIndex)
%     CHANNEL_NUM = 4;
% 
%     for i = 0:CHANNEL_NUM-1
%         FileName = [FileIndex, '_ch_', num2str(i), '.txt'];
%         FileNew = [FilePath, '\', FileName];
%         DataTmp = readtable(FileNew);
%         re_part = double(DataTmp.Var1);
%         im_part = double(DataTmp.Var2);
%         Data(i+1,:) = re_part+1j*im_part;
%     end
% end

function [Data] = SDRDataRead(File)
    DataTmp = readtable(File);
    Ch1 = double(DataTmp.Var1)+1j*double(DataTmp.Var2);
    Ch2 = double(DataTmp.Var3)+1j*double(DataTmp.Var4);
    Ch3 = double(DataTmp.Var5)+1j*double(DataTmp.Var6);
    Ch4 = double(DataTmp.Var7)+1j*double(DataTmp.Var8);
    Data = [Ch1'; Ch2'; Ch3'; Ch4'];
end

