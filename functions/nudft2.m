%{
John Bodenschatz
Marquette University
Rowe Lab
03/22/2023
%}

%{
nudft2.m returns the 2-D nonuniform discrete fourier transform. Assumes
square matrix for image space.

INPUTS:
    x, y (real double): x and y coordinates of image space
    kx, ky (real double): x and y coordinates of k-space
    F (complex double): image space matrix correspinding to x and y inputs
OUTPUT:
    Y (complex vector): fourier transform of F onto kx and ky.
%}

function nudft2 = nudft2(x,y,kx,ky,F)
    x=x(:); y=y(:);
    kx=kx(:); ky=ky(:);
    i = sqrt(-1);
    N = length(kx);
    nudft2 = zeros(N,1);
    F = F(:);
    
    %map x and y onto [0,1], kx and ky onto [0,N]
    x = rescale(x);
    y = rescale(y);
    kx = rescale(kx,0,sqrt(length(F)));
    ky = rescale(ky,0,sqrt(length(F)));

    for j=1:N
        nudft2(j) = sum(F .* exp(-(2*pi*i)*(kx(j)*x + ky(j)*y)));
    end
end