function [] = tsunami(varargin)
%TSUNAMI Summary of this function goes here
%   Parameters: input, output, amplitude, omega, theta, video, videooutput, timestep, maxt, m

input_filename = 'meshes/output.su2';
output_filename = 'data/output.dat';
video_output_filename = 'videos/output.mp4';
amplitude = 1;
omega = 1;                      % [rad/s]
m = 8;                          % Order of truncation of 4.11.1 and 4.11.2 in textbook
theta_I = pi/8;
write_video = false;
time_step = 0.1;
t_end = inf;

if ~isempty(varargin)
    if rem(length(varargin), 2)
        error('tsunami:unevenArgumentCount', 'Error: Uneven argument count. Arguments should follow the "''-key'', value" format. Exiting.');
    end
    for i = 1:2:length(varargin)
        key = varargin{i};
        value = varargin{i+1};

        switch lower(key)
            case "input"
                input_filename = value;
            case "output"
                output_filename = value;
            case "amplitude"
                amplitude = value;
            case "omega"
                omega = value;
            case "theta"
                theta_I = value;
            case "video"
                write_video = value;
            case "videooutput"
                video_output_filename = value;
            case "timestep"
                time_step = value;
            case "maxt"
                t_end = value;
            case "m"
                m = value;
            otherwise
                warning('Warning, unknown parameter: ''%s'', ignoring.', key);
        end
    end
end

if (length(omega) ~= length(amplitude)) || (length(omega) ~= length(theta_I))
    error('tsunami:numberOfWavesNotEqual', 'Error: The number of amplitudes, omega and theta input is not the same. Exiting.');
end

% File input
[points, elements, wall, farfield] = read_su2(input_filename);

% Plotting the mesh
figure()
% Aspect ratio will always be wrong here, as there is no "axis equal" call for 3D plots.
trimesh(elements', points(1, :)', points(2, :)', -points(3, :)');

g = 9.81;                           % [m/s^2]   Gravity
c_water = 1481;                     % [m/s]     Speed of sound
frequency = omega./(2*pi);          % [1/s]     Frequency of wave
period = 1./frequency;              % [s]       Period of wave
lambda = c_water .* period;         % [m]       Wavelength of wave
[R, I] = max(abs(points(1, :)));    % [m]       Farfield radius
h = points(3, I);                   % [m]       Farfield depth
%k = 2*pi./lambda;                   % [1/m]     Wave number 
k = omega./sqrt(g * h);             %           k used in textbook, see 4.1.11 
N_waves = length(omega);            %           Number of waves input

% Printing parameters
for wave = 1:N_waves
    fprintf('Wave with following characteristics:\n');
    fprintf('    Amplitude = %g m\n', amplitude(wave));
    fprintf('    Omega = %g rad/s\n', omega(wave));
    fprintf('    Frequency = %g Hz\n', frequency(wave));
    fprintf('    Period = %g s\n', period(wave));
    fprintf('    Wavelength = %g m\n', lambda(wave));
    fprintf('    Wave number = %g rad/m\n', k(wave));
    fprintf('    Sound speed = %g m/s\n', c_water);
end

M = 2*m + 1;
P = size(farfield, 2);          % Number of nodes on the boundary
E = size(points, 2);            % Total number of nodes in and on the boundary
N_elements = size(elements, 2); % Number of elements

eta = zeros(E, N_waves);        % Displacement

[data_path, data_filename, data_ext] = fileparts(output_filename);

for wave = 1:N_waves
    K_1 = zeros(E, E); % E x E
    K_2 = zeros(M, M); % M x M
    K_3 = zeros(P, M); % P x M

    Q_4 = zeros(P, 1); % 1 x P      transposed (?)
    Q_5 = zeros(M, 1); % 1 x M      transposed (?)

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
                    K_elem(i, j) = K_elem(i, j) - omega(wave)^2/g * Delta_e/6;
                else
                    K_elem(i, j) = K_elem(i, j) - omega(wave)^2/g * Delta_e/2;
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
    K_2(1, 1) = 2 * besselh_prime(0, k(wave)*R) * besselh(0, k(wave)*R);
    for j = 1:m
        K_2(2*j, 2*j) = besselh_prime(j, k(wave)*R) * besselh(j, k(wave)*R);
        K_2(2*j + 1, 2*j + 1) = besselh_prime(j, k(wave)*R) * besselh(j, k(wave)*R);
    end
    K_2 = K_2 * pi * k(wave) * R * h;

    % Building K_3
    L = zeros(P, 1);
    theta = zeros(P, 1);
    q = zeros(P, 1);
    for i = 1:P
        L(i) = sqrt((points(1, farfield(2, i)) - points(1, farfield(1, i)))^2 + (points(2, farfield(2, i)) - points(2, farfield(1, i)))^2);
        sph = to_sph([(points(1, farfield(1, i)) + points(1, farfield(2, i)))/2, (points(2, farfield(1, i)) + points(2, farfield(2, i)))/2, 0]);
        theta(i) = sph(3);
        q(i) = 1i * cos(theta(i) - theta(1)) * exp(1i * k(wave) * R * cos(theta(i) - theta_I(wave)));
    end

    % First row
    K_3(1, 1) = 2 * besselh_prime(0, k(wave)*R) * L(1);
    for j = 1:m
        K_3(1, 2*j) = besselh_prime(j, cos(j * theta(P) + cos(j * theta(1)))) * L(1);
        K_3(1, 2*j+1) = besselh_prime(j, sin(j * theta(P) + sin(j * theta(1)))) * L(1);
    end
    % Other rows
    for i = 2:P
        K_3(i, 1) = 2 * besselh_prime(0, k(wave)*R) * L(i);
        for j = 1:m
            K_3(i, 2*j) = besselh_prime(j, cos(j * theta(i-1) + cos(j * theta(i)))) * L(i);
            K_3(i, 2*j+1) = besselh_prime(j, sin(j * theta(i-1) + sin(j * theta(i)))) * L(i);
        end
    end
    K_3 = -k(wave) * h/2 * K_3;

    % Building Q_4
    % First row
    Q_4(1) = (q(P) - q(1)) * L(1);
    % Other rows
    for i = 2:P
        Q_4(i) = (q(i-1) + q(i)) * L(i);
    end 
    Q_4 = k(wave) * h/2 * Q_4;

    % Building Q_5
    Q_5(1) = besselj(0, k(wave)*R) * besselh_prime(0, k(wave)*R); % Assumes J_0 is J_0(kR)
    for j = 1:m
        Q_5(2*j) = 1i^j * besselj(j, k(wave)*R) * besselh_prime(j, k(wave)*R) * cos(j * theta_I(wave));
        Q_5(2*j + 1) = 1i^j * besselj(j, k(wave)*R) * besselh_prime(j, k(wave)*R) * sin(j * theta_I(wave));
    end
    Q_5 = amplitude(wave) * 2 * pi * R * k(wave) * h * Q_5; %%% CHECK amplitude not here in textbook, but seems like it should be based on 4.11.1 and Q_5 definition

    %% Full system
    LHS = K_3 * (K_2^-1) * (K_3');
    LHS_exp = zeros(size(points, 2), size(points, 2));
    for i = 1:P
        for j = 1:P
            LHS_exp(farfield(1, i), farfield(1, j)) = LHS(i, j);
        end
    end

    K = K_1 - LHS_exp;
    B = Q_4 + K_3*K_2^-1 * Q_5;

    B_exp = zeros(size(points, 2), 1);
    for i = 1:P
        B_exp(farfield(1, i)) = B(i);
    end

    % K * eta = B
    eta(:, wave) = K\B_exp;

    if any(isnan(eta(:, wave)))
        eta(isnan(eta(:, wave)), wave) = 0.0;
        warning('tsunami:nanInEta', 'Warning: There were NaNs in the calculated eta. Usually found in the boundary. Setting to 0.');
    end

    write_solution([data_path filesep data_filename sprintf('_omega%g', omega(wave)) data_ext], eta(:, wave), amplitude(wave), omega(wave));
end

if write_video
    writerObj = VideoWriter(video_output_filename, 'MPEG-4');
    writerObj.FrameRate = 60;
    open(writerObj);
end

if write_video
    figure('Position', [10 10 1920 1080]);
else
    figure();
end

t = 0;
while t <= t_end
    t = t + time_step;
    ksi = zeros(E, 1);
    for wave = 1:N_waves
        ksi = ksi + real(eta(:, wave) * exp(-1i * omega(wave) * t));
    end
    trimesh(elements', points(1, :)', points(2, :)', -points(3, :)');
    hold on
    trimesh(elements', points(1, :)', points(2, :)', ksi);
    hold off
    title(sprintf('%gs', t));
    axis([-R, R, -R, R, -h, 2 * max(amplitude)]);
    drawnow;
    if write_video
        writeVideo(writerObj, getframe(gcf));
    end
end
if write_video
    close(writerObj);
end
end

function [result] = besselh_prime(nu, z)
    result = nu * besselh(nu, z)/z - besselh(nu+1, z);
end