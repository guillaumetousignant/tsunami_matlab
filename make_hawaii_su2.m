function [] = make_hawaii_su2(varargin)
%MAKE_HAWAII Summary of this function goes here
%   Detailed explanation goes here

filename_in = 'data/image.png';
filename_out = 'output.su2';
N_points_farfield = 128;

if ~isempty(varargin)
    if rem(length(varargin), 2)
        error('make_circle_su2:unevenArgumentCount', 'Error, uneven argument count. Arguments should follow the "''-key'', value" format. Exiting.');
    end
    for i = 1:2:length(varargin)
        key = varargin{i};
        value = varargin{i+1};

        switch lower(key)
            case "input"
                filename_in = value;
            case "output"
                filename_out = value;
            case "nptsff"
                N_points_farfield = value;
            otherwise
                warning('Warning, unknown parameter: ''%s'', ignoring.', key);
        end
    end
end

img = imread(filename_in);

figure()
imshow(img);

size_x = size(img, 2);
size_y = size(img, 1);

img_hsv = rgb2hsv(img);

fprintf('Image size: %d, %d\n', size_x, size_y);

radius = min(size_x, size_y)/2;

hold on
center_x = size_x/2 + (size_x/2 - radius);
center_y = size_y/2;
%circle(center_x, center_y, radius);

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

fprintf('Click on points in increasing order on the colour chart and input depth as a positive number. Press any other button when done.\n');
[x,y, button] = ginput(1);
while button == 1
    depth = input('Input depth:\n');
    depth_map_value = [depth_map_value; depth];
    depth_map_hue = [depth_map_hue; img_hsv(round(y), round(x), 1)];
    [x,y, button] = ginput(1);
end

fprintf('\nDepth map:\n');
fprintf('    Depth    Hue\n');
for i = 1:length(depth_map_value)
    fprintf('    %5.5g    %5.5g\n', depth_map_value(i), depth_map_hue(i));
end

points = zeros(0, 3);
%fprintf('Click on a point to get depth.\n');
%[x,y] = ginput(1);
%hue = img_hsv(round(y), round(x), 1);
%disp(hue);
%depth = interp1(depth_map_hue, depth_map_value, hue);
%disp(depth);

ang = linspace(0, 2 * pi, N_points_farfield); 
x_boundary = radius_m * cos(ang);
y_boundary = radius_m * sin(ang);

[x_boundary_pix, y_boundary_pix] = m_to_pixels(x_boundary, y_boundary);

plot(x_boundary_pix, y_boundary_pix);

end

