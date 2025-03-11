clear, clc, close all;

%% Particle Trajectories Plots

% Check if the file exists and is not empty
particleDataFile = 'build/particle_data.txt';
if exist(particleDataFile, 'file') == 2 && dir(particleDataFile).bytes > 0
    % Read the data from the file
    data = readtable(particleDataFile, 'Delimiter', ' ', 'ReadVariableNames', false);

    % Assign column names for easier access
    data.Properties.VariableNames = {'Time', 'ParticleID', 'X', 'Y', 'Z', 'U', 'V', 'W'};

    % Get unique particle IDs
    particleIDs = unique(data.ParticleID);

    % Initialize a structure to hold the data for each particle
    particles = struct();

    % Organize data by particle ID
    for i = 1:length(particleIDs)
        pid = particleIDs(i);
        particleData = data(data.ParticleID == pid, :);
        particles(i).ID = pid;
        particles(i).Time = particleData.Time;
        particles(i).X = particleData.X;
        particles(i).Y = particleData.Y;
        particles(i).Z = particleData.Z;
        particles(i).U = particleData.U;
        particles(i).V = particleData.V;
        particles(i).W = particleData.W;
    end

    % Calculate minimum separation (only for 2 particles)
    if length(particleIDs) == 2
        % Calculate separation at each timestep
        dx = particles(1).X - particles(2).X;
        dy = particles(1).Y - particles(2).Y;
        dz = particles(1).Z - particles(2).Z;
        
        % Calculate distance
        distances = sqrt(dx.^2 + dy.^2 + dz.^2);
        
        % Find minimum separation and its time
        [min_separation, min_idx] = min(distances);
        time_of_min = particles(1).Time(min_idx);
        
        % Display results
        fprintf('Minimum separation between particles: %.6f\n', min_separation);
        fprintf('Time of minimum separation: %.6f\n', time_of_min);
    end

    % Create a 3D plot for the trajectories
    figure;
    hold on;
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Particle Trajectories');

    % Plot the trajectory of each particle
    for i = 1:length(particles)
        plot3(particles(i).X, particles(i).Y, particles(i).Z, '-o', 'DisplayName', ['Particle ' num2str(particles(i).ID)]);
    end

    % Add a legend
    legend('show');
    hold off;
else
    disp('The particle_data.txt file is empty or does not exist. No particle trajectories to plot.');
end

%% Kinetic Energy Plots

% Check if the kinetic energy file exists and is not empty
kineticEnergyFile = 'build/kinetic_energy.txt';
if exist(kineticEnergyFile, 'file') == 2
    % Read the kinetic energy data from the file
    data = readtable(kineticEnergyFile, 'Delimiter', ' ', 'ReadVariableNames', false);

    % Assign column names for easier access
    data.Properties.VariableNames = {'Time', 'KineticEnergy'};

    % Plot the kinetic energy over time
    figure;
    plot(data.Time, data.KineticEnergy, '-o');
    grid on;
    xlabel('Time');
    ylabel('Kinetic Energy');
    title('Kinetic Energy Over Time');
else
    disp('The kinetic_energy.txt file is empty or does not exist. No kinetic energy data to plot.');
end