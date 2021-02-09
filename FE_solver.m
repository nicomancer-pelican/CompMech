function FE_solver(hObject, eventdata, handels, varagin)
    % setup window with grid, inputs panel and axes for results plot
    fig = uifigure('Name', 'AERO96015 - FE Solver');
    fig.Position = [200 200 1000 500];
    grid1 = uigridlayout(fig,[2 1]);
    grid1.RowHeight = {'1x',350};
    p = uipanel(grid1,'Title','FE Solver');
    ax = uiaxes(grid1);
    ax.XGrid = 'on'; ax.YGrid = 'on';

    % grid in the inputs panel
    grid2 = uigridlayout(p,[4 3]);
    grid2.RowHeight = {22,22,22,22,22,22};
    grid2.ColumnWidth = {'1x','1x','1x','1x','1x','1x'};

    % solver label
    solver = uilabel(grid2);
    solver.HorizontalAlignment = 'right';
    solver.Text = 'Solver';

    % solver drop-down
    solver = uidropdown(grid2);
    solver.Items = {'Linear', 'Geometrically Non-Linear'};

    % beam length label
    L = uilabel(grid2);
    L.HorizontalAlignment = 'right';
    L.Text = 'Beam length (m)';

    % beam length field
    L = uieditfield(grid2, 'numeric');
    L.Value = 1;

    % axial stiffness label
    EA = uilabel(grid2);
    EA.HorizontalAlignment = 'right';
    EA.Text = 'Axial stiffness, EA (N)';

    % axial stiffness field
    EI = uieditfield(grid2, 'numeric');
    EI.Value = 10e5;

    % bending stiffness label
    EI = uilabel(grid2);
    EI.HorizontalAlignment = 'right';
    EI.Text = 'Bending stiffness, EI (N m^2)';

    % bending stiffness field
    EI = uieditfield(grid2, 'numeric');
    EI.Value = 10e2;

    % distributed force label
    q = uilabel(grid2);
    q.HorizontalAlignment = 'right';
    q.Text = 'Distributed force (kN/m)';

    % distributed force field
    q = uieditfield(grid2, 'numeric');
    q.Value = 1;

    % no. elements label
    N = uilabel(grid2);
    N.HorizontalAlignment = 'right';
    N.Text = 'Number of elements';

    % no. elements field
    N = uieditfield(grid2, 'numeric');
    N.Value = 1;

    % Create a push button
    btn = uibutton(fig,'push',...
                   'ButtonPushedFcn', @(btn,event) linearFE(btn, L.Value, EA.Value, EI.Value, q.Value, N.value));
end
