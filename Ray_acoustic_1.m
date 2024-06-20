%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Parameters %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z_src = -20;    % Source depth
surface_z = 0;  % Surface location
H = -100;       % Depth
W = 1000;       % Maximum range
theta_aperture = pi/4;  % Aperture angle (in radians)
num_rays = 100;  % Number of rays to plot

source = [0, Z_src]; % Source location
absorption = 0.1;    % Absorption parameter (e.g., dB per meter)
% !!! Ajouter absorption fonction de la frequence

% Updated sound speed profile
depths   = [   0,  -10,  -20,  -35,  -45];
c_values = [1500, 1490, 1480, 1450, 1420];
c = @(z) interp1(depths, c_values, z, 'linear', 'extrap'); % Interpolated sound speed function

figure()
plot(c_values, depths)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Ray %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize intensity range for color mapping
log_intensity_min = inf;
log_intensity_max = -inf;

% Plot rays and compute intensities
rays_r = cell(1, num_rays);
rays_z = cell(1, num_rays);
rays_intensity = cell(1, num_rays);
theta_values = linspace(-theta_aperture/2, theta_aperture/2, num_rays);

for i = 1:num_rays          % For each rays
    theta = theta_values(i); % Initial angle at the source
    [r, z, intensity] = ray_trace(theta, Z_src, H, c, W, absorption);
    rays_r{i} = r;
    rays_z{i} = z;
    rays_intensity{i} = intensity;
    log_intensity = log10(intensity);
    log_intensity_min = min(log_intensity_min, min(log_intensity));
    log_intensity_max = max(log_intensity_max, max(log_intensity));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Plot %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot setup
figure;
hold on;
xlabel('Range (m)');
ylabel('Depth (m)');
title('2D Ray Acoustic Modeling with Depth-Dependent Sound Speed');
axis([0 W H surface_z]);

% Plot surface and bottom
plot([0, W], [surface_z, surface_z], 'k', 'LineWidth', 2);  % Surface
plot([0, W], [H, H], 'k', 'LineWidth', 2);  % Bottom

% Color mapping based on log intensity
colormap('jet');
for i = 1:num_rays
    r = rays_r{i};
    z = rays_z{i};
    intensity = rays_intensity{i};
    log_intensity = log10(intensity);
    for j = 1:length(r)-1
        color_idx = floor((log_intensity(j) - log_intensity_min) / (log_intensity_max - log_intensity_min) * 255) + 1;
        color = jet(256);
        color = color(color_idx, :);  % Get color from colormap based on log intensity
        line([r(j), r(j+1)], [z(j), z(j+1)], 'Color', color);  % Apply color to line
    end
end

% Add color bar for intensity
caxis([log_intensity_min, log_intensity_max]);
colorbar;
ylabel(colorbar, 'Log Intensity');

hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%% Function %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ray tracing function with intensity and absorption
function [r, z, intensity] = ray_trace(theta, Z_src, H, c, W, absorption)
    % Initialize the ray path and intensity
    r = 0;
    z = Z_src;
    intensity = 1;  % Initial intensity
    dt = 0.01;      % Time step
    
    while r(end) < W  % While the ray is in the retangle
        % Calculate sound speed at current depth
        current_c = c(z(end));
        
        % Update ray direction using Snell's law
        if length(r) > 1
            dz = z(end) - z(end-1);
            dr = r(end) - r(end-1);
            theta = atan2(dz, dr);
            next_c = c(z(end) + dz);
            theta = asin((current_c / next_c) * sin(theta));
        end
        
        % Calculate increments in r and z
        dr = current_c * cos(theta) * dt;
        dz = current_c * sin(theta) * dt;
        
        % Update positions
        r = [r, r(end) + dr];
        z = [z, z(end) + dz];
        
        % Update intensity based on absorption
        distance = sqrt(dr^2 + dz^2);
        intensity = [intensity, intensity(end) * 10^(-absorption * distance / 10)];
        
        % Check for reflection on surface
        if z(end) > 0
            z(end) = 2 * 0 - z(end);
            theta = -theta;
        end
        
        % Check for reflection on bottom
        if z(end) < H
            z(end) = 2 * H - z(end);
            theta = -theta;
        end
    end
end
