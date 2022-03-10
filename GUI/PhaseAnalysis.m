function [PhaseDiff] = PhaseAnalysis(DataArray)
%     PhaseDiffArray=DataArray(2:end,:).*conj(DataArray(1,:));
%     [MaxPhaseDiff] = max(angle(PhaseDiffArray),[], 2);
    PhaseDiff = angle(DataArray(1,:)*DataArray(2:end,:)');
end

