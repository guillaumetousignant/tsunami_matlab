close all
clearvars;
clc;

filename_in = 'data/hawaii4.png';
filename_out = 'meshes/hawaii.su2';
N_points_ff = 128;
depth_wall = 32;
ff_ext_factor = 0.975; % Farfield extrude factor
wall_ext_factor = 1.05;
N_domain_r = 32;
N_domain_theta = N_points_ff;
domain_r_exponent = 0.75;
saturation_cutoff = 0.1;
r_start = 0.01;
ratio_high = 2;

img = imread(filename_in);

figure()
imshow(img);

size_x = size(img, 2);
size_y = size(img, 1);

img_hsv = rgb2hsv(img);

fprintf('Image size: %d, %d\n', size_x, size_y);

radius = min(size_x, size_y)/2 - 16;

hold on
center_x = size_x/2 + (size_x/2 - radius) - 16;
center_y = size_y/2;

%% Distance
fprintf('Click on two points with known distance on map.\n');
[x,y] = ginput(2);
distance = input('Input distance in meters:\n');
m2pix = sqrt((x(2) - x(1))^2 + (y(2) - y(1))^2)/distance;
pix2m = 1/m2pix;

pixels_to_m = @(x, y) (deal((x - center_x)*pix2m, (y - center_y)*pix2m));
m_to_pixels = @(x, y) (deal(x * m2pix + center_x, y * m2pix + center_y));

radius_m = radius * pix2m;
disp(radius_m);

depth_map_value = zeros(0, 1);
depth_map_hue = zeros(0, 1);

%% Depth map
fprintf('Click on points in increasing order on the colour chart and input depth as a positive number. Press any other button when done.\n');
[x,y, button] = ginput(1);
while button == 1
    depth = input('Input depth:\n');
    depth_map_value = [depth_map_value; depth];
    depth_map_hue = [depth_map_hue; img_hsv(ceil(y), ceil(x), 1)];
    [x,y, button] = ginput(1);
end

fprintf('\nDepth map:\n');
fprintf('    Depth    Hue\n');
for i = 1:length(depth_map_value)
    fprintf('    %5.5g    %5.5g\n', depth_map_value(i), depth_map_hue(i));
end

%fprintf('Click on a point to get depth.\n');
%[x,y] = ginput(1);
%hue = img_hsv(round(y), round(x), 1);
%disp(hue);
%depth = interp1(depth_map_hue, depth_map_value, hue);
%disp(depth);

%% Farfield
ang = linspace(0, 2 * pi, N_points_ff); 
x_boundary = radius_m * cos(ang);
y_boundary = radius_m * sin(ang);

[x_boundary_pix, y_boundary_pix] = m_to_pixels(x_boundary, y_boundary);

plot([x_boundary_pix, x_boundary_pix(1)], [y_boundary_pix, y_boundary_pix(1)], 'Linewidth', 2, 'Color', 'r', 'LineStyle', '-', 'Marker', 'o');

depths_ff = zeros(N_points_ff, 1);
for i = 1:N_points_ff
    depths_ff(i) = interp1(depth_map_hue, depth_map_value, img_hsv(ceil(y_boundary_pix(i)), ceil(x_boundary_pix(i)), 1), 'linear', 'extrap'); %%% CHECK can nan
end

depths_ff = depths_ff(~isnan(depths_ff)); % Is we fall on black hue will be NaN.
depth_ff = mean(depths_ff);

points_ff = [x_boundary', y_boundary', ones(N_points_ff, 1) * depth_ff];

elements_ff = zeros(N_points_ff, 2);
elements_ff(:, 1) = 1:N_points_ff;
elements_ff(end, 2) = 1;
elements_ff(1:end-1, 2) = 2:N_points_ff;

%% Walls
N_islands = input('Input number of islands:\n');

points_walls = cell(N_islands, 1);
elements_walls = cell(N_islands, 1);
wall_offset = zeros(N_islands+1, 1);
wall_offset(1) = N_points_ff;
N_points_walls = zeros(N_islands, 1);
center_walls = zeros(N_islands, 2);
for i = 1:N_islands
    fprintf('Click on points around island #%d. Press enter when done.\n', i);
    [x,y, button] = ginput();

    plot([x; x(1)], [y; y(1)], 'Linewidth', 2, 'Color', 'b', 'LineStyle', '-', 'Marker', 'o');
    N_points_walls(i) = length(x);
    [x_m, y_m] = pixels_to_m(x, y);
    points_walls{i, 1} = [x_m, y_m, ones(N_points_walls(i), 1)*depth_wall];

    elements_walls{i, 1} = zeros(N_points_walls(i), 1);
    elements_walls{i, 1}(:, 1) = 1+wall_offset(i):N_points_walls(i)+wall_offset(i);
    elements_walls{i, 1}(end, 2) = 1+wall_offset(i);
    elements_walls{i, 1}(1:end-1, 2) = 2+wall_offset(i):N_points_walls(i)+wall_offset(i);

    wall_offset(i + 1) = wall_offset(i) + N_points_walls(i);
    center_walls(i, 1) = mean(x_m);
    center_walls(i, 2) = mean(y_m);
    [center_x, center_y] = m_to_pixels(center_walls(i, 1), center_walls(i, 2));
    plot(center_x, center_y, '+', 'Linewidth', 2, 'Color', 'b');
end

%% Farfield extrude
N_points_ff_ext = N_points_ff;
points_ff_ext = zeros(N_points_ff_ext, 3);

for i = 1:N_points_ff_ext
    point1 = points_ff(elements_ff(i, 1), :);
    point2 = points_ff(elements_ff(i, 2), :);
    x = ff_ext_factor * (point1(1) + point2(1))/2;
    y = ff_ext_factor * (point1(2) + point2(2))/2;
    [x_pix, y_pix] = m_to_pixels(x, y);
    z = (depth_ff + interp1(depth_map_hue, depth_map_value, img_hsv(ceil(y_pix), ceil(x_pix), 1), 'linear', 'extrap'))/2; %%% CHECK can nan
    points_ff_ext(i, :) = [x, y, z];
end

[x_ff_ext_pix, y_ff_ext_pix] = m_to_pixels(points_ff_ext(:, 1), points_ff_ext(:, 2));
plot([x_ff_ext_pix; x_ff_ext_pix(1)], [y_ff_ext_pix; y_ff_ext_pix(1)], 'o', 'Linewidth', 2, 'Color', 'r');

radius_m_ext = ff_ext_factor * radius_m;

%% Islands extrude
N_points_walls_ext = zeros(N_islands, 1);
points_walls_ext = cell(N_islands, 1);

for i = 1:N_islands
    N_points_walls_ext(i, 1) = N_points_walls(i) * 2;
    points_walls_ext{i, 1} = zeros(N_points_walls_ext(i, 1), 3);

    for j = 1:N_points_walls(i)
        point1 = points_walls{i, 1}(elements_walls{i, 1}(j, 1) - wall_offset(i), :);
        point2 = points_walls{i, 1}(elements_walls{i, 1}(j, 2) - wall_offset(i), :);

        x = (point1(1) + point2(1))/2;
        x = center_walls(i, 1) + (x - center_walls(i, 1)) * wall_ext_factor;
        y = (point1(2) + point2(2))/2;
        y = center_walls(i, 2) + (y - center_walls(i, 2)) * wall_ext_factor;
        [x_pix, y_pix] = m_to_pixels(x, y);
        z = interp1(depth_map_hue, depth_map_value, img_hsv(ceil(y_pix), ceil(x_pix), 1), 'linear', 'extrap'); %%% CHECK can nan
        if img_hsv(ceil(y_pix), ceil(x_pix), 2) < saturation_cutoff
            z = NaN;
        end
        points_walls_ext{i, 1}(2 * (j - 1) + 1, :) = [x, y, z];

        x2 = center_walls(i, 1) + (point1(1) - center_walls(i, 1)) * wall_ext_factor^2;
        y2 = center_walls(i, 2) + (point1(2) - center_walls(i, 2)) * wall_ext_factor^2;
        [x2_pix, y2_pix] = m_to_pixels(x2, y2);
        z2 = interp1(depth_map_hue, depth_map_value, img_hsv(ceil(y2_pix), ceil(x2_pix), 1), 'linear', 'extrap'); %%% CHECK can nan
        if img_hsv(ceil(y2_pix), ceil(x2_pix), 2) < saturation_cutoff
            z2 = NaN;
        end
        points_walls_ext{i, 1}(2 * (j - 1) + 2, :) = [x2, y2, z2];
    end

    points_walls_ext{i, 1} = points_walls_ext{i, 1}(~isnan(points_walls_ext{i, 1}(:, 3)), :);
    N_points_walls_ext(i, 1) = size(points_walls_ext{i, 1}, 1);

    [x_wall_ext_pix, y_wall_ext_pix] = m_to_pixels(points_walls_ext{i, 1}(:, 1), points_walls_ext{i, 1}(:, 2));
    plot([x_wall_ext_pix; x_wall_ext_pix(1)], [y_wall_ext_pix; y_wall_ext_pix(1)], 'o', 'Linewidth', 2, 'Color', 'b');
end

%% All other points
points_domain = zeros(N_domain_r * N_domain_theta, 3);
r = linspace(r_start, 1, N_domain_r).^domain_r_exponent * radius_m_ext * ff_ext_factor; % To make them equally distributed
theta = linspace(0, 2 * pi, N_domain_theta);

for j = 1:N_domain_r
    for i = 1:N_domain_theta
        index = i + (j - 1) * N_domain_theta;
        point = to_xyz([r(j), pi/2, theta(i)]);
        [x_pix, y_pix] = m_to_pixels(point(1), point(2));
        point(3) = interp1(depth_map_hue, depth_map_value, img_hsv(ceil(y_pix), ceil(x_pix), 1), 'linear', 'extrap'); %%% CHECK can nan
        if img_hsv(ceil(y_pix), ceil(x_pix), 2) < saturation_cutoff
            point(3) = NaN;
        end
        points_domain(index, :) = point;
    end
end

points_domain = points_domain(~isnan(points_domain(:, 3)), :);
N_points_domain = size(points_domain, 1);

[x_domain_pix, y_domain_pix] = m_to_pixels(points_domain(:, 1), points_domain(:, 2));
plot([x_domain_pix; x_domain_pix(1)], [y_domain_pix; y_domain_pix(1)], 'x', 'Linewidth', 2, 'Color', 'g');

%% All points together now
N_points_wall = sum(N_points_walls);
points_wall = zeros(N_points_wall, 3);
points_wall_ext = zeros(sum(N_points_walls_ext), 3);
offset_ext = 0;
for i = 1:N_islands
    offset = wall_offset(i) - wall_offset(1);
    points_wall(offset + 1:offset + N_points_walls(i), :) = points_walls{i, 1};
    points_wall_ext(1 + offset_ext:N_points_walls_ext(i) + offset_ext, :) = points_walls_ext{i, 1};
    offset_ext= offset_ext + N_points_walls_ext(i);
end

center_point = [0, 0, 0];
[x_center_pix, y_center_pix] = m_to_pixels(center_point(1), center_point(2));
center_point(3) = interp1(depth_map_hue, depth_map_value, img_hsv(ceil(y_center_pix), ceil(x_center_pix), 1), 'linear', 'extrap'); %%% CHECK can nan

points = [points_ff; points_wall; points_ff_ext; points_wall_ext; points_domain; center_point];
points(:, 1) = -points(:, 1); % Images index the other way round oops
N_points = size(points, 1);

%% Delaunay
triangles = delaunay(points(:, 1), points(:, 2));
N_triangles = size(triangles, 1);

%% Deleting 'land' elements
good_triangles = true(N_triangles, 1);
wall_start = N_points_ff + 1;
wall_end = N_points_ff + N_points_wall;

for i = 1:N_triangles
    good_triangles(i, 1) = ~all((triangles(i, :) >= wall_start) & (triangles(i, :) <= wall_end)); %%% CHECK Could also delete elements which have all points on the same island, or elements whose center is land
end

triangles = triangles(good_triangles, :);
N_triangles = size(triangles, 1);

%% Limiters
for i = 1:N_triangles
    ratio1 = abs(points(triangles(i, 1), 3)/points(triangles(i, 2), 3));
    ratio2 = abs(points(triangles(i, 1), 3)/points(triangles(i, 3), 3));
    ratio3 = abs(points(triangles(i, 2), 3)/points(triangles(i, 3), 3));
    ratio4 = abs(points(triangles(i, 2), 3)/points(triangles(i, 1), 3));
    ratio5 = abs(points(triangles(i, 3), 3)/points(triangles(i, 1), 3));
    ratio6 = abs(points(triangles(i, 3), 3)/points(triangles(i, 2), 3));

    if (ratio1 > ratio_high) && (ratio2 > ratio_high)
        points(triangles(i, 1), 3) = (points(triangles(i, 1), 3) + points(triangles(i, 2), 3) + points(triangles(i, 3), 3))/3;
    elseif (ratio3 > ratio_high) && (ratio4 > ratio_high)
        points(triangles(i, 2), 3) = (points(triangles(i, 1), 3) + points(triangles(i, 2), 3) + points(triangles(i, 3), 3))/3;
    elseif (ratio5 > ratio_high) && (ratio6 > ratio_high)
        points(triangles(i, 3), 3) = (points(triangles(i, 1), 3) + points(triangles(i, 2), 3) + points(triangles(i, 3), 3))/3;
    end
end

for i = 1:N_points
    points(i, 3) = max(points(i, 3), depth_wall);
end

%% Plotting
figure()
trisurf(triangles, points(:, 1), points(:, 2), -points(:, 3));

%% Writing file
su2_file = fopen(filename_out, 'w');
    
fprintf(su2_file, 'NDIME= 3\n\n');
fprintf(su2_file, 'NPOIN= %d\n', N_points);
for k = 1:N_points
    fprintf(su2_file, '%g %g %g\n', points(k, 1), points(k, 2), points(k, 3));
end

fprintf(su2_file, '\nNELEM= %d\n', N_triangles);
for k = 1:N_triangles
    fprintf(su2_file, '5 %d %d %d\n', triangles(k, 1), triangles(k, 2), triangles(k, 3));
end

fprintf(su2_file, 'NMARK= %d\n', N_islands + 1);

fprintf(su2_file, 'MARKER_TAG= farfield\n');
fprintf(su2_file, 'MARKER_ELEMS= %d\n', N_points_ff);
for k = 1:N_points_ff
    fprintf(su2_file, '3 %d %d\n',elements_ff(k, 1), elements_ff(k, 2));
end

for i = 1:N_islands
    fprintf(su2_file, 'MARKER_TAG= wall\n');
    fprintf(su2_file, 'MARKER_ELEMS= %d\n', N_points_walls(i));
    for k = 1:N_points_walls(i)
        fprintf(su2_file, '3 %d %d\n', elements_walls{i, 1}(k, 1), elements_walls{i, 1}(k, 2));
    end
end

fclose(su2_file);