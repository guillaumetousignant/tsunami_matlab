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

g = 9.81;                           % [m/s^2] Gravity
c = 1481;                           % [m/s] Speed of sound (NOT SURE)
lambda = c/omega;                   % [m] Wavelength
k = 2*pi/lambda;                    % [1/m] Wave number 
[R, I] = max(abs(points(1, :)));    % [m] Farfield radius
h = points(3, I);

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

% Building K_1
for e = 1:N_elements
    % Element variables
    K_elem = zeros(3, 3);
    x = [points(1, elements(1, e));
         points(1, elements(2, e));
         points(1, elements(3, e))];
    y = [points(2, elements(1, e));
         points(2, elements(2, e));
         points(2, elements(3, e))];
    %a = [x(2)*y(3) - x(3)*y(2); 
    %     x(3)*y(1) - x(1)*y(3);
    %     x(1)*y(2) - x(2)*y(1)];
    b = [y(2) - y(3);
         y(3) - y(1);
         y(1) - y(2)];
    c = [x(3) - x(2);
         x(1) - x(3);
         x(2) - x(1)];
    h_e = points(3, elements(1, e)) + ...
          points(3, elements(2, e)) + ...
          points(3, elements(3, e));

    Delta_e = 0.5 * det([1, x(1), y(1); 1, x(2), y(2); 1, x(3), y(3)]); % Area of element e

    % Assembling element stiffness matrix
    for i = 1:3
        for j = 1:3
            K_elem(i, j) = h_e/(12 * Delta_e) * (b(i)*b(j) + c(i)*c(j));
            if i == j
                K_elem(i, j) = K_elem(i, j) - omega^2/g * Delta_e/6;
            else
                K_elem(i, j) = K_elem(i, j) - omega^2/g * Delta_e/2;
            end
        end
    end

    % Assembling global stiffness matrix
    for i = 1:3
        for j = 1:3
            K_1(elements(i, e), elements(j, e)) = K_1(elements(i, e), elements(j, e)) + K_elem(i, j);
        end
    end
end

% Building K_2
K_2(1, 1) = 2 * besselh_prime(0, k*R) * besselh(0, k*R);
for e = 1:m
    K_2(2*m, 2*m) = besselh_prime(e, k*R) * besselh(e, k*R);
    K_2(2*m + 1, 2*m + 1) = besselh_prime(e, k*R) * besselh(e, k*R);
end
K_2 = K_2 * pi * k * R * h;

% Building K_3
for i = 1:P
    L = sqrt((points(1, farfield(2, i)) - points(1, farfield(1, i)))^2 + (points(2, farfield(2, i)) - points(2, farfield(1, i)))^2);
    sph_1 = to_sph([points(1, farfield(1, i)), points(2, farfield(1, i)), 0]);
    sph_2 = to_sph([points(1, farfield(2, i)), points(2, farfield(2, i)), 0]);
    theta = [sph_1(3);
             sph_2(3)];

    K_3(i, 1) = 2 * besselh_prime(0, k*R) * L;
    for j = 1:m
        K_3(i, 2*m) = besselh_prime(j, cos(j * theta(1) + cos(j * theta(2)))) * L;
        K_3(i, 2*m+1) = besselh_prime(j, sin(j * theta(1) + sin(j * theta(2)))) * L;
    end
end
K_3 = -k * h/2 * K_3;

% Building Q_4
for i = 1:P

end 
Q_4 = k * h/2 * Q_4;

% Building Q_5

K = K_1 - K_3 * (K_2^-1) * (K_3');
B = Q_4 + K_3*K_2^-1 * Q_5;

% K * eta = B

end

function [result] = besselh_prime(nu, z)
    result = nu * besselh(nu, z)/z - besselh(nu+1, z);
end