%{
John Bodenschatz
Marquette University
Rowe Lab
04/18/2023
%}

%{
stats_ttest.m takes in a time series of fmri image data, a task indication 
vector, and performs a t-test on the average task image against the average
non-task image to determine activation

INPUTS:
    imageTS (complex double): 3-D array [Ix, Iy, time] of fmri task
        experiment. assumes coil sensitivities have already been taken care
        of
    design (logical vector): 1s and 0s same length as 3rd dimension of
        kspaceTS. 1s => task.

OUTPUT:
    tmap (real double): [Ix, Iy] map of t-values at each voxel 
%}

function tmap = stats_ttest(imageTS,design)
    nx = size(imageTS,1); ny = size(imageTS,2);
    tmap = zeros(nx, ny);

    %optional masking to remove small t values
    %mask = zeros(nx, ny);
    %mask_threshold = (0.3)*max(max(mean(abs(imageTS),3)));

    % separate into task/nontask images
    restImages=imageTS(:,:,design==0);
    taskImages=imageTS(:,:,design==1);
    restMean = mean(abs(restImages),3);

    figure('name','restmean');
    imagesc(restMean)
    figure('name','taskmean');
    imagesc(mean(abs(taskImages),3))
        
    for jj=1:nx
        for kk=1:ny
            [h,p,ci,stats] = ttest(abs(taskImages(jj,kk,:)),abs(restMean(jj,kk)),'Tail','right');
            tmap(jj,kk) = abs(stats.tstat);
            %mask(jj,kk) = mean(abs(reconstructed_images(jj,kk,:))) > mask_threshold;
        end
    end
end
