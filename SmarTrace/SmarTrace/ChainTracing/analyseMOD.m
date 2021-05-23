function analyseMOD(input,output,LowerLimit,exclude,ko)
SmTr_ChainSamplingMOD(input,exclude)
clc
SmTr_AnalysisMOD(['Re1-randperm-50_10_200_',input],output,LowerLimit,ko)
end