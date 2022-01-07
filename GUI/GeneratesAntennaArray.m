function [Position] = GeneratesAntennaArray(Num,Radius,ArrayType)
    Position = zeros(Num, 3);
    if(strcmp(ArrayType, 'Բ��'))
        RadianUnit = 2*pi/Num;
        RadianTmp = 0:RadianUnit:2*pi;
        for i = 1:Num
            Radian = wrapToPi(RadianTmp(i));
            Position(i,1) = cos(Radian)*Radius;
            Position(i,2) = sin(Radian)*Radius;
        end
    else(strcmp(ArrayType, '����'));
        LenNum = sqrt(Num);
        if(round(LenNum) ~= LenNum)
            message('��������޷�����');
            exit();
        end
    end
end

