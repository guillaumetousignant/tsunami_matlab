function make_circle_su2(varargin)
%MAKE_CIRCLE_SU2 This program creates a mesh with a circular island.
%   The input is made with parameter-value pairs, using the following syntax: 
%   "make_circle_su2('filename', 'meshes/mesh.su2', 'radius', 0.5, 'PARAMETER', VALUE, [...]);"
%   All parameters are optional, defaults are set otherwise.
%
%       Parameter   [default value]         Descrpition 
%
%       filename    ['meshes/input.su2']    Output filename or path to mesh file to be written.
%       nlevels     [8]                     Number of elements in the radial direction.
%       nsides      [32]                    Number of elements in the angular direction.
%       radius      [0.5]                   Radius oh the island, in m.
%       ffradius    [5]                     Radius of the domain, in m.
%       fact        [1.75]                  Bias of the elements towards the island. 1 makes the elements equally spaced radially. Higher than 1 create smaller elements closer to the island, lower than 1 makes them bigger.
%       hfarfield   [10]                    Water depth at the edge of the domain, in m.
%       hwall       [1]                     Water depth at the island, in m.
%       hexp        [0.5]                   Bias of the water depth. 1 is linear, higher than 1 makes it deeper closer to the island, lower than 1 makes it deeper near the edge of the domain.

    % Default values
    filename = 'meshes/input.su2';
    nlevels = 8;
    nsides = 32;
    radius = 0.5;
    ffradius = 5;
    fact = 1.75;
    hfarfield = 10;
    hwall = 1;
    hexp = 0.5;

    % Input parsing
    if ~isempty(varargin)
        if rem(length(varargin), 2)
            error('make_circle_su2:unevenArgumentCount', 'Error, uneven argument count. Arguments should follow the "''-key'', value" format. Exiting.');
        end
        for i = 1:2:length(varargin)
            key = varargin{i};
            value = varargin{i+1};

            switch lower(key)
                case "filename"
                    filename = value;
                case "nlevels"
                    nlevels = value;
                case "nsides"
                    nsides = value;
                case "radius"
                    radius = value;
                case "ffradius"
                    ffradius = value;
                case "fact"
                    fact = value;
                case "hfarfield"
                    hfarfield = value;
                case "hwall"
                    hwall = value;
                case "hexp"
                    hexp = value;
                otherwise
                    warning('Warning, unknown parameter: ''%s'', ignoring.', key);
            end
        end
    end

    npoints = (nlevels + 1) * nsides;
    
    points = zeros(npoints, 3); 
    faces = zeros(nsides * nlevels*2, 3);
    
    angle_step = 2* pi / nsides;
    
    for i = 1:nsides
        point = to_xyz([radius, (i-1)*angle_step, 0]);
        points(i, 1) = point(1);
        points(i, 2) = point(3);
        points(i, 3) = hwall;
    end
    
    for i = 1:nlevels
        for j = 1:nsides
            radius_offset = ((i-1)/(nlevels-1)).^fact * (ffradius - radius) + radius;
            point = to_xyz([radius_offset, (j-1 + 0.5*i)*angle_step, 0]);
            
            points(i*nsides + j, 1) = point(1);
            points(i*nsides + j, 2) = point(3);
            points(i*nsides + j, 3) = ((radius_offset - radius)/(ffradius - radius)).^hexp * (hfarfield - hwall) + hwall;         
            
            temp_point = (i-1)*nsides + j*(j ~= nsides) + 1;
            temp_point2 = i*nsides + j*(j ~= nsides) + 1;
            faces((i-1)*nsides *2 + (j*2)-1, :) = [i*nsides + j, temp_point2, temp_point];

            faces((i-1)*nsides *2 + (j*2), :) = [i*nsides + j, temp_point, (i-1)*nsides + j];
        end        
    end

    plot_circle();
    plot_circle_3D();
    write_su2(filename);

    function plot_circle()
        figure()
        for k = 1:size(faces,1)
            x = [points(faces(k, 1), 1), points(faces(k, 2), 1), points(faces(k, 3), 1), points(faces(k, 1), 1)];
            y = [points(faces(k, 1), 2), points(faces(k, 2), 2), points(faces(k, 3), 2), points(faces(k, 1), 2)];
            plot(x, y);
            hold on
        end
        axis equal
    end

    function plot_circle_3D()
        figure()
        trimesh(faces, points(:, 1), points(:, 2), -points(:, 3));
        xlim([-ffradius, ffradius]);
        ylim([-ffradius, ffradius]);
        zlim([-2*ffradius, 0]);
    end

    function write_su2(filename)
        su2_file = fopen(filename, 'w');
    
        fprintf(su2_file, 'NDIME= 3\n\n');
        fprintf(su2_file, 'NPOIN= %g\n', npoints);
        for k = 1:npoints
            fprintf(su2_file, '%g %g %g\n', points(k, 1), points(k, 2), points(k, 3));
        end
    
        fprintf(su2_file, '\nNELEM= %g\n', nsides * nlevels*2);
        for k = 1:nsides * nlevels*2
            fprintf(su2_file, '5 %g %g %g\n', faces(k, 1), faces(k, 2), faces(k, 3));
        end

        fprintf(su2_file, 'NMARK= 2\n');
        fprintf(su2_file, 'MARKER_TAG= wall\n');
        fprintf(su2_file, 'MARKER_ELEMS= %g\n', nsides);
        for k = 1:nsides
            fprintf(su2_file, '3 %g %g\n', k, k*(k~=nsides) + 1);
        end

        fprintf(su2_file, 'MARKER_TAG= farfield\n');
        fprintf(su2_file, 'MARKER_ELEMS= %g\n', nsides);
        fprintf(su2_file, '3 %g %g\n', nsides * nlevels + nsides, nsides * nlevels + 1);
        for k = 2:nsides
            fprintf(su2_file, '3 %g %g\n', nsides * nlevels + k - 1, nsides * nlevels + k);
        end

        fclose(su2_file);
    end
end

