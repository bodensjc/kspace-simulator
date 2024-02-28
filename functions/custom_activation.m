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

function actMap = custom_activation(ImX, ImY)
    actMap = zeros(ImY,ImX);

    % Number of squares along X and Y:
    gridSizeX = ImX;
    gridSizeY = ImY;

    % Initialise figure:
    fig = figure('Position', [100 100 500 500], 'Units', 'Pixels');

    % Set squares boundaries inside figure in pixels:
    gridPosition = [0 0 1 1] .* repmat(fig.Position(3:4), 1, 2);

    % Compute square size:
    buttonSizeX = gridPosition(3) / gridSizeX ;
    buttonSizeY = gridPosition(4) / gridSizeY;

    % Get squares coordiantes in pixels:
    [X,Y] = meshgrid(gridPosition(1) : buttonSizeX : sum(gridPosition([1,3])) - buttonSizeX, ...
                     gridPosition(2) : buttonSizeY : sum(gridPosition([2,4])) - buttonSizeY);

    % Loop over the grid just defined:
    for col = 1 : gridSizeX
        for row = 1 : gridSizeY
            % Get position:
            position = [X(row, col), Y(row, col), buttonSizeX, buttonSizeY];
            % Create square as pushbutton:
            uicontrol(fig, ...
                      'BackgroundColor','g', ...
                      'Position', position, ...
                      'Units', 'Pixels', ...
                      'UserData', [row, col, 1], ...
                      'Callback', @click);
        end
    end

    function click(obj,~)
    % Check if button was on or off and change color:
    if obj.UserData(3)
        obj.BackgroundColor = 'r';
        obj.UserData(3) = 0;
    else
        obj.BackgroundColor = 'g';
        obj.UserData(3) = 1;
    end
    end


end