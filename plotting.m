clear,clc, close all;

% Read the data from the file
data = readtable('particle_data.txt', 'Delimiter', ' ', 'ReadVariableNames', false);

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

% Plot the position components over time for each particle
figure;
subplot(3, 1, 1);
hold on;
title('Position X over Time');
xlabel('Time');
ylabel('X');
for i = 1:length(particles)
    plot(particles(i).Time, particles(i).X, '-o', 'DisplayName', ['Particle ' num2str(particles(i).ID)]);
end
legend('show');
hold off;

subplot(3, 1, 2);
hold on;
title('Position Y over Time');
xlabel('Time');
ylabel('Y');
for i = 1:length(particles)
    plot(particles(i).Time, particles(i).Y, '-o', 'DisplayName', ['Particle ' num2str(particles(i).ID)]);
end
legend('show');
hold off;

subplot(3, 1, 3);
hold on;
title('Position Z over Time');
xlabel('Time');
ylabel('Z');
for i = 1:length(particles)
    plot(particles(i).Time, particles(i).Z, '-o', 'DisplayName', ['Particle ' num2str(particles(i).ID)]);
end
legend('show');
hold off;

% Plot the velocity components over time for each particle
figure;
subplot(3, 1, 1);
hold on;
title('Velocity U over Time');
xlabel('Time');
ylabel('U');
for i = 1:length(particles)
    plot(particles(i).Time, particles(i).U, '-o', 'DisplayName', ['Particle ' num2str(particles(i).ID)]);
end
legend('show');
hold off;

subplot(3, 1, 2);
hold on;
title('Velocity V over Time');
xlabel('Time');
ylabel('V');
for i = 1:length(particles)
    plot(particles(i).Time, particles(i).V, '-o', 'DisplayName', ['Particle ' num2str(particles(i).ID)]);
end
legend('show');
hold off;

subplot(3, 1, 3);
hold on;
title('Velocity W over Time');
xlabel('Time');
ylabel('W');
for i = 1:length(particles)
    plot(particles(i).Time, particles(i).W, '-o', 'DisplayName', ['Particle ' num2str(particles(i).ID)]);
end
legend('show');
hold off;