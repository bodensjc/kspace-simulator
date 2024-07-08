%{
John Bodenschatz
Marquette University
Rowe Lab
07/03/2023
%}

%{
kspace_Cartesian.m takes in the relevant MRI information and returns an 
array of k-space positions and the relative times at which they are measured"


INPUTS:
    MRI (struct): MRI struct from SHAKER
    rec_size (int): length of square reconstructed image size

OUTPUT:
    kspaceTrajectory (real double): three column matrix. [kx, ky, time]
%}

function kspaceTrajectory = kspace_Cartesian(MRI,rec_size)
    accelerationFactor = MRI.AccelerationFactor;
    echoTime = MRI.EchoTime;
    

    [kx,ky] = meshgrid(-rec_size/2:rec_size/2-1,-rec_size/2:rec_size/2-1);
    kx = kx(1:accelerationFactor:end,:);
    ky = ky(1:accelerationFactor:end,:);
    %nkx = width(kx); nky = height(kx);


    % create time map
    eesp = 0.000720;
    deltaT = 0;%1/(2*125*10^3);
    extras=0;
    TE=echoTime*1000;
    
    timeMap = zeros(rec_size,rec_size);
    count = 0;
    for row = 1:rec_size
        for col=1:rec_size
            if mod(row,2)==1
                timeMap(row,col)=(row-1)*extras*deltaT+(col-1)*deltaT+(row-1)*eesp;
            elseif mod(row,2)==0
                timeMap(row,rec_size-col+1)=(row-1)*extras*deltaT+(col-1)*deltaT+(row-1)*eesp;
            end
        end
    end
    timeMapMsec = timeMap*1000;
    init_t = TE - timeMapMsec(rec_size/2-1,rec_size);
    timeMap = timeMap + init_t/1000;


end