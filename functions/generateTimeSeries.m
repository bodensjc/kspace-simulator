%{
John Bodenschatz
Marquette University
Rowe Lab
2/29/24
%}

%{
generateTimeSeries.m takes in a kspace with and without activation,
relevant task  design info (HRF, epochs, etc.) and noise stats to return a
complete timem series data ksapce.

INPUTS:
    kspaceR (complex double): rest state kspace image
    kspaceT (complex double): task state kspace image
    nEpochs (integer): number of epochs
    nRest (integer): number of rest state immages per epoch
    nTask (integer): number of task state images per epoch
    sigma (double): st dev of real and imaginary noise for fourier space

OUTPUT:
    kspaceTS (complex double): time series of kspace arrays
%}

function kspaceTS = generateTimeSeries(kspaceR,kspaceT,nEpochs,nRest,nTask,sigma)
    restBlock=repmat(kspaceR,1,1,nRest);
    taskBlock=repmat(kspaceT,1,1,nTask);
    singleBlock=cat(3,restBlock,taskBlock);
    kspaceTS=repmat(singleBlock,1,1,nEpochs);
    kspaceTS = kspaceTS + sigma*randn(size(kspaceTS)) + sigma*randn(size(kspaceTS))*sqrt(-1);
end