function [MaxPhaseDiff, MeanPhaseDiff] = PhaseAnalysis(DataArray)
    PhaseDiffArray=DataArray(2:end,:).*conj(DataArray(1,:));
    [MaxPhaseDiff, Index] = max(angle(PhaseDiffArray),[], 2);
    MeanPhaseDiff = angle((DataArray(1,:)*DataArray(2:end,:)')/size(DataArray, 2));
end

