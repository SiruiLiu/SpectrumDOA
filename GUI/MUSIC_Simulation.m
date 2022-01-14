function [signal_num] = MUSIC_Simulation(varargin)
    SampledData = varargin{1};
    Position = varargin{2};
    Lambda = varargin{3};
    element_num = length(SampledData(:,1));
    snapshot_num = length(SampledData(1,:));
    Covariance = SampledData*SampledData'./snapshot_num;
    angle_step = pi/180;       % ½Ç¶ÈËÑË÷²½½ø 
    theta = 0:angle_step:pi/2;     % ¸©Ñö½ÇËÑË÷·¶Î§
    phi = 0:angle_step:2*pi;     % ·½Î»½ÇËÑË÷·¶Î§
            
%     CovarianceNew = Toeplitz(SampledData);
%     CovarianceNew2 = MMusic(Covariance);
    if(nargin <= 3)
        [O, R] = eig(Covariance);    
        Q = sort(diag(R), 'descend');
        Pmdl = zeros(1, element_num - 1);
        for k = 1:element_num-1
            temp = element_num - k;
            Qk = Q(k + 1 : element_num);
            L(k) = snapshot_num * temp * log(1 / temp * sum(Qk) / power(prod(Qk), 1 / temp));
            Pmdl(k) = k * (2 * element_num - k) * log (snapshot_num / 2);
        end

        Kmdl = L + Pmdl;
        signal_num = find(Kmdl == min(Kmdl));
    else
        signal_num = varargin{3};
    end
    [O, R] = eig(Covariance); 
    O = O(1:element_num, 1:element_num);
    R = diag(R(1:element_num, 1:element_num));

    % Æ×·åËÑË÷
    UU=O(:,1:element_num-(signal_num));
    Pmusic = zeros(length(theta),length(phi));
    for i = 1:length(theta)
        for j = 1:length(phi)
            unit = [sin(theta(i))*cos(phi(j)),sin(theta(i))*sin(phi(j)),cos(theta(i))];
            SteeringVector = exp(1i*2*pi/Lambda*Position*unit');
            WW = SteeringVector'*(UU*UU')*SteeringVector;
            Pmusic(i,j) = abs(1/WW);
        end
    end

    % Pmusic=10*log10(Pmusic/max(max(Pmusic)));

    figure(102);
    meshc(Pmusic);
    grid on
    
    figure(103)   % ÈýÎ¬ÊÓÍ¼
    meshc(sin(theta)'*cos(phi),sin(theta)'*sin(phi),Pmusic);
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

