%{
John Bodenschatz
Marquette University
Rowe Lab
06/26/2023
%}

%{
getnoise.m Takes a given image and signal to noise ratio (SNR) then returns
           the scale parameter of a Rayleigh distribution that can be used
           to generate noise in k-space. The scale parameter, sigma, of the
           Rayleigh dist. is the standard deviation of the related
           bivariate normal distribution.

INPUTS:
    IMG (real/complex double): Image-space representation of object 
    SNR (real double): desired SNR

OUTPUT:
    getnoise (real double): Rayleigh dist. scale parameter sigma
%}

function globalSNR = globalSNR(IMG,MASK)
    IMG = abs(IMG);
    %mask=IMG>0;
    %mask = (IMG./max(IMG,[],'all')) > (1/10); % 1/10 kind of arbitrary...
    % get mean of signal region, sum of nonzero pixels / count
    IMGmean = sum(MASK.*IMG,'all') ./ sum(MASK~=0,'all');
    % get stdev of non signal region
    SPACEmean = sum((1-MASK).*IMG,'all') ./ (numel(IMG)-sum(MASK~=0,'all'));
    SPACEstdev = sqrt(sum(((1-MASK).*IMG-SPACEmean).^2,'all') ./ (numel(IMG)-sum(MASK~=0,'all')));
    % k-space stats, scale image space by dimensions
    globalSNR = IMGmean / SPACEmean;
end
