function [] = tsunami(varargin)
%TSUNAMI Summary of this function goes here
%   Detailed explanation goes here

filename = 'output.su2';
amplitude = 1;
omega = 1;
theta_I = pi/8;

if ~isempty(varargin)
    if rem(length(varargin), 2)
        error('tsunami:unevenArgumentCount', 'Error, uneven argument count. Arguments should follow the "''-key'', value" format. Exiting.');
    end
    for i = 1:2:length(varargin)
        key = varargin{i};
        value = varargin{i+1};

        switch lower(key)
            case "-filename"
                filename = value;
            case "-amplitude"
                amplitude = value;
            case "-omega"
                omega = value;
            case "-theta"
                theta_I = value;
            otherwise
                warning('Warning, unknown parameter: ''%s'', ignoring.', key);
        end
    end
end

% File input
[points, elements, wall, farfield] = read_su2(filename);

% Plotting the mesh
figure()
% Aspect ratio will always be wrong here, as there is no "axis equal" call for 3D plots.
trimesh(elements', points(1, :)', points(2, :)', -points(3, :)');

g = 9.81;                       % [m/s^2] Gravity
c = 1481;                       % [m/s] Speed of sound (NOT SURE)
lambda = c/omega;               % [m] Wavelength
k = 2*pi/lambda;                % [1/m] Wave number 
m = 8;                          % Order of truncation of 4.11.1 and 4.11.2 in textbook
M = 2*m + 1;
P = size(farfield, 2);          % Number of nodes on the boundary
E = size(points, 2);            % Total number of nodes in and on the boundary
N_elements = size(elements, 2); % Number of elements

K_1 = zeros(E, E); % E x E
K_2 = zeros(M, M); % M x M
K_3 = zeros(P, M); % P x M

Q_4 = zeros(1, P); % 1 x P
Q_5 = zeros(1, M); % 1 x M

for k = 1:N_elements
    K_elem = zeros(3, 3);
    for i = 1:3
        for j = 1:3
            K_elem(i, j) = 
        end
    end
end

K = K_1 - K_3 * (K_2^-1) * (K_3');

end

