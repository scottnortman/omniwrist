function plotDebug3d( items, ax_scale, lims )
    %   FUNCTION plotDebug3d( items )
    %       Raises a new figure and constructs an axis suitable for 3d viewing
    %       and debug

    % expets a cell array as an input; 4x4 arrays are assumed to be homogeneous
    % transforms

    fighan = figure( 'name', 'Plot Debug 3D Pro',...
                    'numbertitle', 'off', ...
                    'menubar', 'none' );

    axhan = axes( 'parent', fighan,...
                    'projection', 'orthographic', ...
                    'units', 'normalized', ....
                    'position', [0 0 1 1] );
    axhan.Toolbar.Visible = 'on';
    set( axhan, 'nextplot', 'replacechildren' );
    axis vis3d equal 
    grid on
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    xlim( lims( 1,: ) );
    ylim( lims( 2, : ) );
    zlim( lims( 3, : ) );
    


    for ii=1:length( items )

        hgt = hgtransform( 'parent', axhan );
        
        item = items{ii};

        [nr,nc] = size( item );

        % handle 4x4 as hgt, 4x1 or 3x1 as a point
        if( nr == 4 ) && ( nc == 4 )
            
            plot3( hgt,...
                [0 ax_scale], [0 0], [0 0], 'r',...
                [0 0], [0 ax_scale], [0 0], 'g',...
                [0 0], [0 0], [0 ax_scale], 'b',...
                0, 0, 0, 'k*', 'linewidth', 2 );
            
            hgt.Matrix = item;

        elseif nc==1
            plot3( hgt,...
                0, 0, 0, 'k*', ...
                'linewidth', 2 );
            hgt.Matrix(:,4) = item;
        else
            warning('Skipped index %d', ii );
        end


    end

end