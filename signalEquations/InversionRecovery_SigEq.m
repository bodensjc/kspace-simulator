%{
John Bodenschatz
Marquette University
Rowe Lab
12/05/2023
%}

%{
InversionRecovery_SigEq.m returns the simulated kspace given relevant MRI info
    and kspace sampling points. Using Equation (9)* from the paper
    "Incorporating Relaxivities..."

    INVERSION RECOVERY

INPUTS:
    M0 (real double): proton spin density map
    T2star (real double): T2star map
    T1 (real double): T2star map
    deltaB (real double): magnetic field inhomogeneity map
    timeMap (real double): echo time TE map; can be single number (TE)
    TR (real double): Repitition time; TR~1
    kx, ky (real doubles): x and y coordinates of k-space

OUTPUT:
    kspace (complex double): simulated kspace from M0 and other details
%}

function kspace = InversionRecovery_SigEq(M0,T1,deltaB,kx,ky,MRI)
    gamma = MRI.gamma; % 42 MHz/T (H nuclei gyromagnetic ratio)
    TR=MRI.RepititionTime;
    i=sqrt(-1); % imaginary unit

    % vectorize kspace
    kspace_size = size(kx);
    kxx=kx(:); kyy=ky(:);
    numKpts = length(kxx);
   
    
    img_length=length(M0);
    [x, y] = meshgrid(linspace(0,1,img_length),linspace(0,1,img_length));

    % build the integrand of signal equation
    ksp = @(j) M0.*(1-2*exp(-TI./T1)+exp(-TR./T1)).*...
        exp(-i*2*pi*(kxx(j)*x+kyy(j)*y));
    if deltaB ~= -1
        ksp =@(j) ksp(j).*exp(i*gamma*deltaB); % time map??? 3/1/2024
    end

    kspace=zeros(numKpts,1);
    for j=1:numKpts
        kspace(j) = sum(ksp(j),'all');
    end
    kspace = reshape(kspace,kspace_size);
end