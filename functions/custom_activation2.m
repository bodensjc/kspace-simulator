%{
John Bodenschatz
Marquette University
Rowe Lab
10/29/2023
%}

%{
custom_activation.m includes the code to create a custom activation map for
the purposes of k-space simulation.

INPUTS:
    ImX, ImY (ints): image x and y directions

OUTPUT:
    actMap (real double): matrix of dim [ImY, ImX] that as ones where
    activation is desired, zeroes elsewhere
%}

function actMap = custom_activation2(ImX, ImY)
    %get an example image
    figure,IM=repelem(get(image,'CData'),ImX,ImY);close

    %create figure with image
    h.f=figure('Units','Normalized');
    h.ax=axes('Parent',h.f);
    h.im=imshow(IM,[],'Parent',h.ax);
    %set callbacks
    h.ax.ButtonDownFcn=@Callback;
    h.im.ButtonDownFcn=@Callback;
    %store in guidata
    guidata(h.f,h)
    function Callback(hObject,eventdata)
        %load the guidata handles struct
        h=guidata(hObject);
        clc,get(h.f,'CurrentPoint')
    end


end