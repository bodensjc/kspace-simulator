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

function tmap = stats_ttest(imageTS,design)
    nx = size(imageTS,1); ny = size(imageTS,2);
    SNRmap = zeros(nx, ny);

    %optional masking to remove small t values
    %mask = zeros(nx, ny);
    %mask_threshold = (0.3)*max(max(mean(abs(imageTS),3)));

    % separate into task/nontask images
    restImages=imageTS(:,:,:,design==0);
    taskImages=imageTS(:,:,:,design==1);
    restMean = mean(mean(abs(restImages),4),3);
        
    for jj=1:nx
        for kk=1:ny
            [h,p,ci,stats] = ttest(abs(taskImages(jj,kk,:)),abs(restMean(jj,kk)),'Tail','right');
            SNRmap(jj,kk) = abs(stats.tstat);
            %mask(jj,kk) = mean(abs(reconstructed_images(jj,kk,:))) > mask_threshold;
        end
    end
end
