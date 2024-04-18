addpath('C:\Users\johnb\Documents\fmri-research\Data')
addpath('C:\Users\johnb\Documents\fmri-research\HAPI\\HAPIvsIRGN_CP_ext\Simulations\As')
addpath('C:\Users\johnb\Documents\fmri-research\kspace-simulator\data')
addpath('C:\Users\johnb\Documents\fmri-research\kspace-simulator\functions')
%% load maps and choose slice
load M0Map96.mat
load T1Map96.mat
load T2Map96.mat
slice = 47;
recon_size = length(M0Map);
M0=M0Map(:,:,slice);
T2star=T2Map(:,:,slice);
T1=T1Map(:,:,slice);
deltaB = repmat(linspace(-2.5,2.5,length(T1))*10e-7,length(T1),1);
deltaB = deltaB + T1./(10^7);
i=sqrt(-1);
n_phantom = recon_size;
%% mr parameters
TE = 50e-03;
TR=1;
MR.FOV = 24;% cm 
MR.gym = 4258.57; % Hz/G
MR.BW = 125000; % Hz
MR.G = MR.BW/(MR.FOV * MR.gym);
flipAngle = 90;
gamma = 42.58e06;

%% HAPI settings
accelfactor=1;
HAPI_samples = 128;
n_spokes = 128;%ceil(pi/2*HAPI_samples / accelfactor);
pa_HAPI = linspace(0,180,n_spokes+1); pa_HAPI(end)=[];
npa = length(pa_HAPI);

pa = linspace(-90,90,n_spokes+1);
pa(end)=[];
%% make kspcace
% gridding
kx=zeros(HAPI_samples,npa);
ky=zeros(HAPI_samples,npa);
xr=[-HAPI_samples/2:1:HAPI_samples/2-1]';
for m=1:npa
    kx(:,m) = xr*cosd(pa_HAPI(m));
    ky(:,m) = xr*sind(pa_HAPI(m));
end

[kx,ky]=meshgrid(-64:63,-64:63);

rampdcf = 1/npa + (kx.^2 + ky.^2).^(0.5);

trueImage = M0.*(1-exp(-TR./T1)).*exp(-TE./T2star).*...
        (sind(flipAngle)./(1-cosd(flipAngle)*exp(-TR./T1))).*...
        exp(i*gamma*deltaB*TE);
%k_space = kspace(M0,T1,T2star,deltaB,TE,TR,flipAngle,kx,ky);
%% save images
figure;
imagesc(angle(trueImage))
axis image; colormap(gray); axis off; colorbar;
%exportgraphics(gcf,['hapi-stuff-2-1-24/trueImage-phase-colorbar.png'],'Resolution',600)
%close all
%% object with noise for task activation
load ActMap96.mat

nRepitions = 50;


CNR = 0.75; pctIncrease=0.02; % cnr and signal increase percent
objmax = max(max(abs(trueImage)));
b1 = pctIncrease*objmax;
sigma = b1/CNR;

trueImage = (trueImage./abs(trueImage)) .* (abs(trueImage)+b1*ActMap);

trueImage = repmat(trueImage,1,1,nRepitions);

%% add noise
trueImage = trueImage + sigma*randn(size(trueImage)) + sigma*randn(size(trueImage))*sqrt(-1);

%% coil setup
cp = importdata('coil_profile.mat');
%cp1 = cp(:,:,1);
%cp=cp1;
ncoils = size(cp,3);
cp_sig_res = zeros(n_phantom,n_phantom,ncoils);
cp_rec_res = zeros(recon_size,recon_size,ncoils);
for c=1:ncoils
    cp_sig_res(:,:,c) = imresize(cp(:,:,c),[n_phantom n_phantom]);
    cp_rec_res(:,:,c) = imresize(cp(:,:,c),[recon_size recon_size]);
end
%% HAPI signal generation
[rawdata, projdataHS, projdataRS, k] = signal_generator(trueImage,cp_sig_res,pa,pa_HAPI,HAPI_samples,recon_size,MR);
projHS = reshape(projdataHS,HAPI_samples*ncoils*npa,1);
projRS = reshape(projdataRS,recon_size*ncoils*npa,1);

%projHS = projHS + 0.05*randn(size(projHS)) + sqrt(-1)*0.05*randn(size(projHS));

%% HAPI A generation
tic
A_totHS = zeros(HAPI_samples*ncoils*npa,recon_size^2);
A = zeros(HAPI_samples,recon_size^2);
for ppp = 1:npa
    A = HAPI_core(pa_HAPI(ppp),HAPI_samples,recon_size,ncoils*length(cp_rec_res));
    disp([num2str(ppp), 'th projection'])
    for k = 1:ncoils
        dum = reshape(cp_rec_res(:,:,k),1,recon_size^2);
        dum2 = repmat(dum,HAPI_samples,1);
        A_totHS((ppp-1)*HAPI_samples*ncoils+(k-1)*HAPI_samples+1:(ppp-1)*HAPI_samples*ncoils+k*HAPI_samples,:) = A.*dum2;
    end
end
toc
clear A; clear dum; clear dum2
%% Projection A generation
A_totRS = zeros(recon_size*ncoils*npa,recon_size^2);
A = zeros(recon_size,recon_size^2);
for i = 1:npa
    A = HAPI_core(pa_HAPI(i),recon_size,recon_size);
    disp([num2str(i), 'th projection'])
    for k = 1:ncoils
        dum = reshape(cp_rec_res(:,:,k),1,recon_size^2);
        dum2 = repmat(dum,recon_size,1);
        A_totRS((i-1)*recon_size*ncoils+(k-1)*recon_size+1:(i-1)*recon_size*ncoils+k*recon_size,:) = A.*dum2;
    end
end

%% HAPI reconstruction
tic
HAPI_reconHS = reshape(lsqr(A_totHS,projHS,[]),[recon_size recon_size]);
%HAPI_reconHS = abs(HAPI_reconHS)/norm(abs(HAPI_reconHS));
figure; imagesc(abs(HAPI_reconHS)); colormap gray; axis image; axis off;
%exportgraphics(gcf,'plot.png','Resolution',600)
toc
%% HAPI reconstruction lsqr normal equations
% DOES NOT WORK WELL DO NOT USE
tic
%HAPI_reconHS2 = reshape((A_totHS'*A_totHS)\(A_totHS'*projHS) ,[recon_size recon_size]);
HAPI_reconHS2 = abs(HAPI_reconHS2)/norm(abs(HAPI_reconHS2));
figure; imagesc(abs(HAPI_reconHS2)); colormap gray; axis image; axis off;
toc
%% Projection reconstruction
HAPI_reconRS = reshape(lsqr(A_totRS,projRS,[],30),[recon_size recon_size]);
HAPI_reconRS = abs(HAPI_reconRS)/norm(abs(HAPI_reconRS));
figure; imagesc(abs(HAPI_reconRS)); colormap gray;
%exportgraphics(gcf,'plot.png','Resolution',600)


%% repeated stack
reconstructed_images = zeros(size(trueImage));
for ll=1:nRepitions
    obj = trueImage(:,:,ll);
    [rawdata, projdataHS, projdataRS, k] = signal_generator(obj,cp_sig_res,pa,pa_HAPI,HAPI_samples,recon_size,MR);
    projHS = reshape(projdataHS,HAPI_samples*ncoils*npa,1);
    disp([num2str(ll), 'th image out of ', num2str(nRepitions)])
    tic
    HAPI_reconHS = reshape(lsqr(A_totHS,projHS,[]),[recon_size recon_size]);
    reconstructed_images(:,:,ll) = HAPI_reconHS;
    toc
end

figure;
imagesc(abs(trueImage(:,:,15)))
axis image; colorbar; colormap(gray)    

%% VERY SLOW (1.25 hrs for 25, 151spoke, 96 samples/spoke, 96x96 images to reconstruct)
tic
HAPI_reconHS = A_totHS\projHS;
toc
HAPI_reconHS2 = reshape(HAPI_reconHS,recon_size,recon_size,nRepitions);
figure;
imagesc(abs(HAPI_reconHS2(:,:,1)))
axis image; colorbar; colormap(gray)


%% t testing for activation

activation = zeros(recon_size,recon_size);

reconstructed_images = task;

mask = zeros(recon_size,recon_size);
mask_threshold = (0.3)*max(max(mean(abs(reconstructed_images),3)));

%test = abs(reconstructed_images) - abs(HAPI_reconHS);

for jj=1:recon_size
    for kk=1:recon_size
        [h,p,ci,stats] = ttest(abs(reconstructed_images(jj,kk,:)),abs(notaskMean(jj,kk)),'Tail','right');
        activation(jj,kk) = stats.tstat;
        mask(jj,kk) = mean(abs(reconstructed_images(jj,kk,:))) > mask_threshold;
        %activation(jj,kk) = mean(abs(reconstructed_images(jj,kk,:)));
    end
end
figure;
imagesc(abs(activation))
colorbar; axis image; colormap(hot)

figure;
imagesc(mask)
axis image; colormap(hot)

figure;
imagesc(activation.*mask)
colorbar; axis image; colormap(hot); caxis([0 max(max(activation.*mask))])

%% SNR

variance = zeros(recon_size,recon_size);
SNR = zeros(recon_size,recon_size);
reconstructed_images = task;

for jj=1:recon_size
    for kk=1:recon_size
        SNR(jj,kk) = mean(abs(reconstructed_images(jj,kk,:)))/std(abs(reconstructed_images(jj,kk,:)));
        variance(jj,kk) = var(abs(reconstructed_images(jj,kk,:)));
    end
end
figure;
imagesc(abs(variance))
colorbar; axis image; colormap(hot)
%% save images
figure;
imagesc(activation)
axis image; colormap(hot); axis off; colorbar;
caxis([0 max(max(activation.*mask))])
exportgraphics(gcf,['hapi-stuff-2-1-24/76sp-raw-activation.png'],'Resolution',600)
%close all

%% CODE TO RUN OVERNIGHT
clear all
addpath('C:\Users\johnb\Documents\fmri-research\Data')
addpath('C:\Users\johnb\Documents\fmri-research\HAPI\\HAPIvsIRGN_CP_ext\Simulations\As')
addpath('C:\Users\johnb\Documents\fmri-research\kspace-simulator\data')
addpath('C:\Users\johnb\Documents\fmri-research\kspace-simulator\functions')

load A_totHS_96samp_151sp_equal.mat

load M0Map96.mat
load T1Map96.mat
load T2Map96.mat
load ActMap96.mat
slice = 47;
recon_size = length(M0Map);
M0=M0Map(:,:,slice);
T2star=T2Map(:,:,slice);
T1=T1Map(:,:,slice);
deltaB = repmat(linspace(-2.5,2.5,length(T1))*10e-7,length(T1),1);
deltaB = deltaB + T1./(10^8);
i=sqrt(-1);
n_phantom = recon_size;

accelfactor=1;
HAPI_samples = 96;
n_spokes = ceil(pi/2*HAPI_samples / accelfactor);
pa_HAPI = linspace(0,180,n_spokes+1); pa_HAPI(end)=[];
npa = length(pa_HAPI);
pa = linspace(-90,90,n_spokes+1);
pa(end)=[];

TE = 50e-03;
TR=1;
MR.FOV = 24;% cm 
MR.gym = 4258.57; % Hz/G
MR.BW = 125000; % Hz
MR.G = MR.BW/(MR.FOV * MR.gym);
flipAngle = 90;
gamma = 42.58e06;
cp = importdata('coil_profile.mat');
ncoils = size(cp,3);
cp_sig_res = zeros(n_phantom,n_phantom,ncoils);
cp_rec_res = zeros(recon_size,recon_size,ncoils);
for c=1:ncoils
    cp_sig_res(:,:,c) = imresize(cp(:,:,c),[n_phantom n_phantom]);
    cp_rec_res(:,:,c) = imresize(cp(:,:,c),[recon_size recon_size]);
end

%% no task
disp('START 96SAMP 151SPOKE 50 IMAGES NO TASK')
trueImage = M0.*(1-exp(-TR./T1)).*exp(-TE./T2star).*...
        (sind(flipAngle)./(1-cosd(flipAngle)*exp(-TR./T1))).*...
        exp(i*gamma*deltaB*TE);
nRepitions = 50;
CNR = 0.75; pctIncrease=0.02; % cnr and signal increase percent
objmax = max(max(abs(trueImage)));
b1 = pctIncrease*objmax;
sigma = b1/CNR;
%trueImage = (trueImage./abs(trueImage)) .* (abs(trueImage)+b1*ActMap);
trueImage = repmat(trueImage,1,1,nRepitions);
trueImage = trueImage + sigma*randn(size(trueImage)) + sigma*randn(size(trueImage))*sqrt(-1);

reconstructed_images = zeros(size(trueImage));
for ll=1:nRepitions
    obj = trueImage(:,:,ll);
    [rawdata, projdataHS, projdataRS, k] = signal_generator(obj,cp_sig_res,pa,pa_HAPI,HAPI_samples,recon_size,MR);
    projHS = reshape(projdataHS,HAPI_samples*ncoils*npa,1);
    disp([num2str(ll), 'th image out of ', num2str(nRepitions)])
    tic
    HAPI_reconHS = reshape(lsqr(A_totHS,projHS,[]),[recon_size recon_size]);
    reconstructed_images(:,:,ll) = HAPI_reconHS;
    toc
end
save('96samp_151spoke_recon_noise_notask.mat','reconstructed_images','-v7.3')
disp('DONE WITH 96SAMP 151SPOKE 50 IMAGES NO TASK')


% task
disp('START 96SAMP 151SPOKE 50 IMAGES TASK')
trueImage = M0.*(1-exp(-TR./T1)).*exp(-TE./T2star).*...
        (sind(flipAngle)./(1-cosd(flipAngle)*exp(-TR./T1))).*...
        exp(i*gamma*deltaB*TE);
nRepitions = 50;
CNR = 0.75; pctIncrease=0.02; % cnr and signal increase percent
objmax = max(max(abs(trueImage)));
b1 = pctIncrease*objmax;
sigma = b1/CNR;
trueImage = (trueImage./abs(trueImage)) .* (abs(trueImage)+b1*ActMap);
trueImage = repmat(trueImage,1,1,nRepitions);
trueImage = trueImage + sigma*randn(size(trueImage)) + sigma*randn(size(trueImage))*sqrt(-1);

reconstructed_images = zeros(size(trueImage));
for ll=1:nRepitions
    obj = trueImage(:,:,ll);
    [rawdata, projdataHS, projdataRS, k] = signal_generator(obj,cp_sig_res,pa,pa_HAPI,HAPI_samples,recon_size,MR);
    projHS = reshape(projdataHS,HAPI_samples*ncoils*npa,1);
    disp([num2str(ll), 'th image out of ', num2str(nRepitions)])
    tic
    HAPI_reconHS = reshape(lsqr(A_totHS,projHS,[]),[recon_size recon_size]);
    reconstructed_images(:,:,ll) = HAPI_reconHS;
    toc
end
save('96samp_151spoke_recon_noise_task.mat','reconstructed_images','-v7.3')
disp('DONE WITH 96SAMP 151SPOKE 50 IMAGES TASK')


clear all







%% REPEAT FOR 16 SPOKES, 512 samples



clear all
addpath('C:\Users\johnb\Documents\fmri-research\Data')
addpath('C:\Users\johnb\Documents\fmri-research\HAPI\\HAPIvsIRGN_CP_ext\Simulations\As')
addpath('C:\Users\johnb\Documents\fmri-research\kspace-simulator\data')
addpath('C:\Users\johnb\Documents\fmri-research\kspace-simulator\functions')

load A_totHS_512samp_16sp_equal.mat

load M0Map96.mat
load T1Map96.mat
load T2Map96.mat
load ActMap96.mat
slice = 47;
recon_size = length(M0Map);
M0=M0Map(:,:,slice);
T2star=T2Map(:,:,slice);
T1=T1Map(:,:,slice);
deltaB = repmat(linspace(-2.5,2.5,length(T1))*10e-7,length(T1),1);
deltaB = deltaB + T1./(10^8);
i=sqrt(-1);
n_phantom = recon_size;

accelfactor=1;
HAPI_samples = 512;
n_spokes = 16;
pa_HAPI = linspace(0,180,n_spokes+1); pa_HAPI(end)=[];
npa = length(pa_HAPI);
pa = linspace(-90,90,n_spokes+1);
pa(end)=[];

TE = 50e-03;
TR=1;
MR.FOV = 24;% cm 
MR.gym = 4258.57; % Hz/G
MR.BW = 125000; % Hz
MR.G = MR.BW/(MR.FOV * MR.gym);
flipAngle = 90;
gamma = 42.58e06;
cp = importdata('coil_profile.mat');
ncoils = size(cp,3);
cp_sig_res = zeros(n_phantom,n_phantom,ncoils);
cp_rec_res = zeros(recon_size,recon_size,ncoils);
for c=1:ncoils
    cp_sig_res(:,:,c) = imresize(cp(:,:,c),[n_phantom n_phantom]);
    cp_rec_res(:,:,c) = imresize(cp(:,:,c),[recon_size recon_size]);
end



% no task
disp('START 512SAMP 16SPOKE 50 IMAGES NO TASK')
trueImage = M0.*(1-exp(-TR./T1)).*exp(-TE./T2star).*...
        (sind(flipAngle)./(1-cosd(flipAngle)*exp(-TR./T1))).*...
        exp(i*gamma*deltaB*TE);
nRepitions = 50;
CNR = 0.75; pctIncrease=0.02; % cnr and signal increase percent
objmax = max(max(abs(trueImage)));
b1 = pctIncrease*objmax;
sigma = b1/CNR;
%trueImage = (trueImage./abs(trueImage)) .* (abs(trueImage)+b1*ActMap);
trueImage = repmat(trueImage,1,1,nRepitions);
trueImage = trueImage + sigma*randn(size(trueImage)) + sigma*randn(size(trueImage))*sqrt(-1);

reconstructed_images = zeros(size(trueImage));
for ll=1:nRepitions
    obj = trueImage(:,:,ll);
    [rawdata, projdataHS, projdataRS, k] = signal_generator(obj,cp_sig_res,pa,pa_HAPI,HAPI_samples,recon_size,MR);
    projHS = reshape(projdataHS,HAPI_samples*ncoils*npa,1);
    disp([num2str(ll), 'th image out of ', num2str(nRepitions)])
    tic
    HAPI_reconHS = reshape(lsqr(A_totHS,projHS,[]),[recon_size recon_size]);
    reconstructed_images(:,:,ll) = HAPI_reconHS;
    toc
end
save('512samp_16spoke_recon_noise_notask.mat','reconstructed_images','-v7.3')
disp('DONE WITH 512SAMP 16SPOKE 50 IMAGES NO TASK')


% task
disp('START 512SAMP 16SPOKE 50 IMAGES TASK')
trueImage = M0.*(1-exp(-TR./T1)).*exp(-TE./T2star).*...
        (sind(flipAngle)./(1-cosd(flipAngle)*exp(-TR./T1))).*...
        exp(i*gamma*deltaB*TE);
nRepitions = 50;
CNR = 0.75; pctIncrease=0.02; % cnr and signal increase percent
objmax = max(max(abs(trueImage)));
b1 = pctIncrease*objmax;
sigma = b1/CNR;
trueImage = (trueImage./abs(trueImage)) .* (abs(trueImage)+b1*ActMap);
trueImage = repmat(trueImage,1,1,nRepitions);
trueImage = trueImage + sigma*randn(size(trueImage)) + sigma*randn(size(trueImage))*sqrt(-1);

reconstructed_images = zeros(size(trueImage));
for ll=1:nRepitions
    obj = trueImage(:,:,ll);
    [rawdata, projdataHS, projdataRS, k] = signal_generator(obj,cp_sig_res,pa,pa_HAPI,HAPI_samples,recon_size,MR);
    projHS = reshape(projdataHS,HAPI_samples*ncoils*npa,1);
    disp([num2str(ll), 'th image out of ', num2str(nRepitions)])
    tic
    HAPI_reconHS = reshape(lsqr(A_totHS,projHS,[]),[recon_size recon_size]);
    reconstructed_images(:,:,ll) = HAPI_reconHS;
    toc
end
save('512samp_16spoke_recon_noise_task.mat','reconstructed_images','-v7.3')
disp('DONE WITH 512SAMP 16SPOKE 50 IMAGES TASK')
