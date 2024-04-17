%{
John Bodenschatz
Marquette University
Rowe Lab
10/19/2023
%}

%{
kspace_GE.m returns the simulated kspace given relevant MRI info
    and kspace sampling points. Using Equation (9)* from the paper
    "Incorporating Relaxivities..."

    GRADIENT ECHO (MAIN FMRI)

INPUTS:
    M0 (real double): proton spin density map
    T2star (real double): T2star map
    T1 (real double): T2star map
    deltaB (real double): magnetic field inhomogeneity map
    timeMap (real double): echo time TE map; can be single number (TE)
    TR (real double): Repitition time; TR~1
    kx, ky (real doubles): x and y coordinates of k-space
    coilSensitivity (real double): coil sensitivity array

OUTPUT:
    kspace (complex double): simulated kspace from M0 and other details
%}

function kspace = kspace_GE(M0,T1,T2star,deltaB,timeMap,TR,flipAngle,kx,ky,coilSensitivity)
    gamma = 42.58e06; % 42 MHz/T (H nuclei gyromagnetic ratio)
    i=sqrt(-1); % imaginary unit

    % vectorize kspace
    nCoils = size(coilSensitivity,3);
    kspace_size = size(kx);
    kspace=zeros(size(kx,1),size(kx,2),nCoils);
    kxx=kx(:); kyy=ky(:);
    numKpts = length(kxx);
    timeMap = timeMap(:);

    if prod(size(timeMap),'all')==1 % if just TE is given instead of a map
        timeMap = repmat(timeMap,numKpts,1);
    end
    
    img_length=length(M0);
    [x, y] = meshgrid(linspace(0,1,img_length),linspace(0,1,img_length));

    
    for c=1:nCoils
        % build the integrand of signal equation
        ksp = @(j) coilSensitivity(:,:,c).*M0.*(1-exp(-TR./T1)).*exp(-timeMap(j)./T2star).*...
            (sind(flipAngle)./(1-cosd(flipAngle)*exp(-TR./T1))).*...
            exp(-i*2*pi*(kxx(j)*x+kyy(j)*y));
        if deltaB ~= -1
            ksp =@(j) ksp(j).*exp(i*gamma*deltaB*timeMap(j));
        end
    
        kspaceC=zeros(numKpts,1);
        for j=1:numKpts
            kspaceC(j) = sum(ksp(j),'all');
        end
        kspaceC = reshape(kspaceC,kspace_size);
        kspace(:,:,c)=kspaceC;
    end
end