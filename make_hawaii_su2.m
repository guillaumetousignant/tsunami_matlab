close all
clearvars;
clc;

filename_in = 'data/hawaii2.png';
filename_out = 'meshes/hawaii.su2';
N_points_ff = 128;
depth_wall = 32;
ff_ext_factor = 0.975; % Farfield extrude factor
wall_ext_factor = 1.1;

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
    depths_ff(i) = interp1(depth_map_hue, depth_map_value, img_hsv(ceil(y_boundary_pix(i)), ceil(x_boundary_pix(i)), 1)); %%% CHECK can nan
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
    z = (depth_ff + interp1(depth_map_hue, depth_map_value, img_hsv(ceil(y_pix), ceil(x_pix), 1)))/2; %%% CHECK can nan
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
        z = interp1(depth_map_hue, depth_map_value, img_hsv(ceil(y_pix), ceil(x_pix), 1)); %%% CHECK can nan
        points_walls_ext{i, 1}(2 * j, :) = [x, y, z];

        x2 = center_walls(i, 1) + (point1(1) - center_walls(i, 1)) * wall_ext_factor^2;
        y2 = center_walls(i, 2) + (point1(2) - center_walls(i, 2)) * wall_ext_factor^2;
        [x2_pix, y2_pix] = m_to_pixels(x2, y2);
        z2 = interp1(depth_map_hue, depth_map_value, img_hsv(ceil(y2_pix), ceil(x2_pix), 1)); %%% CHECK can nan
        points_walls_ext{i, 1}(2 * j + 1, :) = [x2, y2, z2];
    end

    [x_wall_ext_pix, y_wall_ext_pix] = m_to_pixels(points_walls_ext{i, 1}(:, 1), points_walls_ext{i, 1}(:, 2));
    plot([x_wall_ext_pix; x_wall_ext_pix(1)], [y_wall_ext_pix; y_wall_ext_pix(1)], 'o', 'Linewidth', 2, 'Color', 'b');
end