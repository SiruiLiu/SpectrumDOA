function [Position] = GeneratesAntennaArray(Num,Radius,ArrayType)
    Position = zeros(Num, 3);
    if(strcmp(ArrayType, '圆阵'))
        RadianUnit = 2*pi/Num;
        RadianTmp = 0:RadianUnit:2*pi;
        for i = 1:Num
            Radian = wrapToPi(RadianTmp(i));
            Position(i,1) = cos(Radian)*Radius;
            Position(i,2) = sin(Radian)*Radius;
        end
    else(strcmp(ArrayType, '方阵'));
        LenNum = sqrt(Num);
        if(round(LenNum) ~= LenNum)
            message('矩阵个数无法开方');
            exit();
        end
    end
end

