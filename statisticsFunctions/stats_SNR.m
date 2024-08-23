%{
John Bodenschatz
Marquette University
Rowe Lab
04/18/2023
%}

%{
stats_SNR.m takes in a time series of fmri image data, a task indication 
vector, and performs a t-test on the average task image against the average
non-task image to determine activation

INPUTS:
    imageTS (complex double): 4-D array [Ix, Iy, coils, time] of fmri task
        experiment.
    design (logical vector): 1s and 0s same length as 3rd dimension of
        kspaceTS. 1s => task.

OUTPUT:
    tmap (real double): [Ix, Iy] map of t-values at each voxel 
%}

function SNRmap = stats_SNR(imageTS,design,TaskRest)
    nx = size(imageTS,1); ny = size(imageTS,2);

    %optional masking to remove small t values
    %mask = zeros(nx, ny);
    %mask_threshold = (0.3)*max(max(mean(abs(imageTS),3)));

    % separate into task/nontask images
    restImages=squeeze(mean(abs(imageTS(:,:,:,design==0)),3));
    taskImages=squeeze(mean(abs(imageTS(:,:,:,design==1)),3));
    if TaskRest == 0
        restMean = squeeze(mean(restImages,3));
        SNRmap = restMean./sqrt(var(restImages,0,3));
    elseif TaskRest == 1
        taskMean = squeeze(mean(taskImages,3));
        SNRmap = taskMean./sqrt(var(taskImages,0,3));
    end
end
