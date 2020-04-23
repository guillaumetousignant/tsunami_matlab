function [points, elements, wall, farfield] = read_su2(filename)
%READ_SU2 Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(filename,'r');

if fid == -1
    error('read_su2:invalidInput', ['Unable to open file "', filename, '".']);
end

fgetl(fid); % NDIM
fgetl(fid); % Spacing
tline = fgetl(fid); % NPOIN
[ln, rest] = strtok(tline);

if ~strcmp(ln, 'NPOIN=')
    error('read_su2:invalidToken', ['Invalid token "', ln, '". Expected "NPOIN=".']);
end

n_points = str2double(rest);
points = zeros(3, n_points);

for i = 1:n_points
    tline = fgetl(fid);
    points(:, i) = sscanf(tline,'%f');
end

tline = fgetl(fid); % Spacing
tline = fgetl(fid); % NELEM
[ln, rest] = strtok(tline);

if ~strcmp(ln, 'NELEM=')
    error('read_su2:invalidToken', ['Invalid token "', ln, '". Expected "NELEM=".']);
end

n_elements = str2double(rest);
elements = zeros(3, n_elements);

for i = 1:n_elements
    tline = fgetl(fid);
    [ln, rest] = strtok(tline);
    if ~strcmp(ln, '5')
        error('read_su2:invalidElementType', ['Invalid element type "', ln, '". Only works for triangles, "5".']); 
    end
    elements(:, i) = sscanf(rest,'%f');
end

tline = fgetl(fid); % NMARK
[ln, rest] = strtok(tline);
if ~strcmp(ln, 'NMARK=')
    error('read_su2:invalidToken', ['Invalid token "', ln, '". Expected "NMARK=".']);
end

n_markers = str2double(rest);
n_walls = n_markers - 1;
wall_index = 1;
wall = cell(n_walls, 1);
farfield = [];

for i = 1:n_markers
    tline = fgetl(fid); % MARKER_TAG
    [ln, rest] = strtok(tline);
    if ~strcmp(ln, 'MARKER_TAG=')
        error('read_su2:invalidToken', ['Invalid token "', ln, '". Expected "MARKER_TAG=".']);
    end

    switch strtrim(lower(rest))
        case 'wall'
            tline = fgetl(fid); % MARKER_ELEMS
            [ln, rest] = strtok(tline);
            if ~strcmp(ln, 'MARKER_ELEMS=')
                error('read_su2:invalidToken', ['Invalid token "', ln, '". Expected "MARKER_ELEMS=".']);
            end

            n_wall = str2double(rest);
            wall{wall_index, 1} = zeros(2, n_wall); 

            for j = 1:n_wall
                tline = fgetl(fid);
                [ln, rest] = strtok(tline);
                if ~strcmp(ln, '3')
                    error('read_su2:invalidElementType', ['Invalid element type "', ln, '". Only works for lines, "3".']); 
                end
                wall{wall_index, 1}(:, j) = sscanf(rest,'%f');
            end
            wall_index = wall_index + 1;

        case 'farfield'
            tline = fgetl(fid); % MARKER_ELEMS
            [ln, rest] = strtok(tline);
            if ~strcmp(ln, 'MARKER_ELEMS=')
                error('read_su2:invalidToken', ['Invalid token "', ln, '". Expected "MARKER_ELEMS=".']);
            end

            n_farfield = str2double(rest);
            farfield = zeros(2, n_farfield); 

            for j = 1:n_farfield
                tline = fgetl(fid);
                [ln, rest] = strtok(tline);
                if ~strcmp(ln, '3')
                    error('read_su2:invalidElementType', ['Invalid element type "', ln, '". Only works for lines, "3".']); 
                end
                farfield(:, j) = sscanf(rest,'%f');
            end

        otherwise
            error('read_su2:invalidMarkerTag', ['Invalid marker tag "', rest, '". Expected "wall" or "farfield".']);
    end
end

if isempty(wall)
    error('read_su2:wallNotFound', 'Wall marker not found.');
end

if isempty(farfield)
    error('read_su2:farfieldNotFound', 'Far field marker not found.');
end

fclose(fid);
end

