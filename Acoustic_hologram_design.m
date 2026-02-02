% 3D transducer simulation using time reversal to design acoustic holograms 
% Author: Johnny Tang
% Date: October 2025
clear;
clc;
addpath(genpath('k-Wave'))



%% Set background medium properties
source_freq = 1.55e6;                                       % Hz
alpha_power = 1.1;
temp = 17;                                                  % degrees celsius
ppw = 5;                                                    % points per wavelength


c_lens  = 1438.81;                                             % speed of sound in lens material (m/s)
d_lens  = 1009.60;                                             % density of lens material (kg/m^3)
alpha_lens = 2.94/(source_freq*10^-6);                       % absorption coefficient of lens material (alpha=alpha0/f^y) y = 1

c_water = waterSoundSpeed(temp);                            % speed of sound of water (m/s)
d_water = waterDensity(temp);                               % density of water (kg/m^3)
alpha_water = waterAbsorption(source_freq*1e-6,temp)/1.68;  % absorption coefficient of water (dB/cmMHz^y)

vars = {'temp'};                                            % clear any variables not needed again
clear(vars{:})

%% Create the computational grid and load data
lambda = c_water/source_freq;   % wavelength (m)
dx = lambda/ppw;  	            % grid point spacing in the x direction (m)
dy = dx;                        % grid point spacing in the y direction (m)
dz = dx;                        % grid point spacing in the z direction (m)

Lx = 0.0905;                    % grid size (m)
Lz = 0.06;

% set number of grid points 
Nx = round(Lx/dx);
Ny = Nx;
Nz = round(Lz/dz);

% create the kgrid structure
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

vars = {'lambda','dx','dy','dz','ppw','Lx','Lz','Nx','Ny','Nz'}; % clear any variables not needed again
clear(vars{:})

%% Generate a binary mask for the H195 transducer
% define the active element of the curved transducer with a hole in the middle
ROC = round(0.045/kgrid.dx);  % 0.06165                   % radius of curvature of the transducer (grid points)

bowl_pos      = [kgrid.Nx/2, kgrid.Ny/2, 1];        % (grid points)
focus_pos     = [kgrid.Nx/2, kgrid.Ny/2, ROC];      % (grid points)

% dimensions of the transducer
d1        = 0.064;       % 0.0842                         % diameter (m)
diameter1 = (2*floor((d1/kgrid.dx)/2)+1);           % round diameter to nearest odd number of grid points
 
% dimensions of the hole in the transducer
d2        = 0.0424;    % 0.0436                            % diameter (m)
diameter2 = (2*floor((d2/kgrid.dx)/2)+1);           % round diameter to nearest odd number of grid points

% radius of curvature
radius    = (2*floor((ROC)/2)+1);                   % round radius to nearest odd number of grid points

% make a binary mask for the transducer
transducer = makeBowl([kgrid.Nx, kgrid.Ny, kgrid.Nz], bowl_pos, radius, diameter1, focus_pos)...
    -makeBowl([kgrid.Nx, kgrid.Ny, kgrid.Nz], bowl_pos, radius, diameter2, focus_pos);

[~,ind_t] = find(transducer(:,round(kgrid.Nx/2),:)==1);
transducer = transducer(:,:,min(ind_t):max(ind_t)); % crop the transducer matrix 

vars = {'ROC','bowl_pos','focus_pos','d1','diameter1','d2','diameter2','radius'}; % clear any variables not needed again
clear(vars{:})

%% Define the holographic plane as a binary mask
% find the holographic plane (the exit plane of the transducer) in the z direction
holographic_planeZ = max(ind_t)-min(ind_t)+1;

% create a binary mask for the holographic plane 
holographic_plane = zeros(kgrid.Nx,kgrid.Ny,kgrid.Nz);
for m = 1:kgrid.Nx
    for n = 1:kgrid.Ny
        x = max(transducer(m,n,:));
        if x~=0
            holographic_plane(m,n,holographic_planeZ) = 1;
        end
    end
end

vars = {'ind_t','m','n','x','transducer_small'};   % clear any variables not needed again
clear(vars{:})

%% Create mask for hologram behind the holographic plane
hologram_back = zeros(kgrid.Nx,kgrid.Ny,kgrid.Nz);
holo_plane2D = holographic_plane(:,:,holographic_planeZ);
for i = 1:kgrid.Nx
    for j = 1:kgrid.Ny
        if holo_plane2D(i,j)==1
            val = holographic_planeZ;
            while val>0
            hologram_back(i,j,val)=1; % set lens voxels to 1
            val = val - 1;
            end
            val_t = round(mean(find(transducer(i,j,:)==1))); % set voxels behind the transducer to 0 
            while val_t>0
            hologram_back(i,j,val_t)=0;
            val_t = val_t - 1;
            end
        end
    end
end

vars = {'holo_plane2D','i','j','val','val_t'};  % clear any variables not needed again
clear(vars{:})

%% Set medium properties
% define the properties of the propagation medium for a homogeneous medium with the hologram backing
medium.sound_speed                              = ones(kgrid.Nx,kgrid.Ny,kgrid.Nz)*c_water;     % (m/s)
medium.sound_speed(hologram_back==1)            = c_lens;                                       

medium.density                                  = ones(kgrid.Nx,kgrid.Ny,kgrid.Nz)*d_water;     % (kg/m^-3)
medium.density(hologram_back==1)                = d_lens;                                       

medium.alpha_coeff                              = ones(kgrid.Nx,kgrid.Ny,kgrid.Nz)*alpha_water; % (dB/(MHz^y*cm)) 
medium.alpha_coeff(hologram_back==1)            = alpha_lens;                                   

medium.alpha_power                              = alpha_power;

% make a time array for the simulations
makeTime(kgrid,medium.sound_speed,0.2)

vars = {'t','height_cone','hologram_back','c_brain','c_skull','d_brain','d_skull','alpha_water','alpha_brain','alpha_skull','alpha_power'};
clear(vars{:})

%% Transducer simulation 
% set the time array length 
kgrid.t_array = 0:kgrid.dt:(holographic_planeZ*kgrid.dx*2)/c_water;
kgrid.Nt = length(kgrid.t_array);

% define a time varying sinusoidal source
source_mag      = 400e3;          % (Pa)
sampling_freq   = 1/kgrid.dt;     % (Hz)
num_cycles      = 200;
source.p        = source_mag*toneBurst(sampling_freq, source_freq, num_cycles,'Envelope','Rectangular');
source.p_mask   = transducer;     % set the transducer as the emitting element

% define the sensor at the holographic plane
sensor.mask = holographic_plane;

% run the simulation
input_args = {'PlotPML', true,'PMLInside', false, 'SaveToDisk', 'kwave_input_data_transducergrey.h5'};
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

sensor_data = h5read('kwave_output_data_transducergrey.h5', '/p');

% extract phase information at the holographic plane
Nfft = ceil(length(sensor_data(1,:))/2)+1;
f = linspace(0,0.5,Nfft)*sampling_freq;
ind = dsearchn(f',source_freq);
Pfilt_aux = fft(sensor_data,[],2);
P_transducer = Pfilt_aux(:,ind);

%% Target simulation
% set the time array length 
focaldepth = round(0.05/kgrid.dx);   % (grid points)
kgrid.t_array = 0:kgrid.dt:((focaldepth-holographic_planeZ)*kgrid.dx*2)/c_water;
kgrid.Nt = length(kgrid.t_array);

% define the target of the hologram as 2 monoples targetdistance apart
targetdistance = 0.01;              % (m)
source.p_mask = zeros(kgrid.Nx,kgrid.Ny,kgrid.Nz);
source.p_mask(round(kgrid.Nx/2),round(kgrid.Ny/2-(targetdistance/2)/kgrid.dy),focaldepth)=1;
source.p_mask(round(kgrid.Nx/2),round(kgrid.Ny/2+(targetdistance/2)/kgrid.dy),focaldepth)=1;

% ADDITION: place a third dot to form a triangular arrangement
% For an equilateral triangle the vertical offset should be:
vertical_offset = targetdistance * sqrt(3) / 2;  % offset in meters
% Place the third dot below (increasing row index) and centered in column:
source.p_mask(round(kgrid.Nx/2) + round(vertical_offset/kgrid.dx), round(kgrid.Ny/2), focaldepth) = 1;

% define the sensor at the holographic plane
sensor.mask = holographic_plane;

% run the simulation
input_args = {'PlotPML', true,'PMLInside', false, 'SaveToDisk', 'kwave_input_data_targetgreytriple.h5'};
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

sensor_data = h5read('kwave_output_data_targetgreytriple.h5', '/p');

% extract phase information at the holographic plane
Nfft = ceil(length(sensor_data(1,:))/2);
f = linspace(0,0.5,Nfft)*sampling_freq;
ind = dsearchn(f',source_freq);
Pfilt_aux = fft(sensor_data,[],2);
P_target = Pfilt_aux(:,ind);

vars = {'sensor_data','Nfft','f','ind','Pfilt_aux'};    % clear any variables not needed again
clear(vars{:})

%% Phase and height calculations
% find the total phase information
P_total = conj(-1*P_target).*conj(P_transducer);
phaseangle = angle(P_total);    

save ('phasetransducer.mat','P_transducer','-v7.3')
save ('phasetarget.mat','P_target','-v7.3')

% convert phase information to height
Z_norm = (d_lens*c_lens)/(d_water*c_water); % normal impedance 
k_medium = 2*pi*source_freq/c_water;        % wave number of medium
k_lens = 2*pi*source_freq/c_lens;           % wave number of lens

% interpolation
h1 = 0:1e-9:c_lens/source_freq; %0.04 for PU;             % interpolation points 
d  = holographic_planeZ*kgrid.dx;           % distance from transducer to holographic plane (m)
T  = (2*Z_norm*exp(-1i*k_medium*(d-h1)))./(2*Z_norm*cos(k_lens*h1)+1i*(Z_norm^2+1)*sin(k_lens*h1)); % transmitted wave equation
T_arg = angle(T);                        

height = interp1(T_arg,h1,phaseangle,"linear","extrap");                % linear interpolation to find the height values
height_reshape = NaN(kgrid.Nx,kgrid.Ny);
height_reshape(holographic_plane(:,:,holographic_planeZ)==1)=height;    % reshape the height values into the hologrpahic plane shape 

h = height_reshape/kgrid.dx;                % pixel heights in the lens (grid points)

save('h.mat','h','-v7.3')

vars = {'P_transducer_rev','P_target','P_transducer','P_total','phaseangle','Z_norm','d_water','k_medium','k_lens','h1','d','T','T_arg','height','holographic_plane','height_reshape'};
clear(vars{:})

%% Hologram surface generation 
% create surface of holographic lens and transducer
[X,Y] = meshgrid((1:kgrid.Nx)*kgrid.dx,(1:kgrid.Ny)*kgrid.dy);

% create matrix of actual height values of hologram (grid points)
hologram_surf = h+holographic_planeZ;

% create matrix of transducer surface height values (grid points)
transducer_surf = zeros(kgrid.Nx,kgrid.Ny);

% allocate the height values of the hologram to a 2x2 matrix
for i = 1:kgrid.Nx
    for j = 1:kgrid.Ny
        [~,ind_t] = max(transducer(i,j,:));
        transducer_surf(i,j) = ind_t;
    end
end

hologram_surf_stl = ones(kgrid.Nx,kgrid.Ny)*max(transducer_surf(:));
transducer_surf_stl = hologram_surf_stl;
hologram_surf_stl(~isnan(hologram_surf))=hologram_surf(~isnan(hologram_surf));
transducer_surf_stl(transducer_surf~=1)=transducer_surf(transducer_surf~=1);

% export STL files
surf2stl('hologram.stl',X*1000,Y*1000,(hologram_surf_stl*kgrid.dx*1000),'ascii')        % convert and save hologram to an stl file in mm
surf2stl('transducer.stl',X*1000,Y*1000,(transducer_surf_stl*kgrid.dx*1000),'ascii')    % convert and save transducer to an stl file in mm

vars = {'X','Y','h','holographic_planeZ','transducer_surf','ind_t','hologram_surf_stl','transducer_surf_stl'};  % clear any variables not needed again
clear(vars{:})

%% Test the hologram
% create a binary map for the holographic lens in a 3D matrix 
hologram_map = zeros(kgrid.Nx,kgrid.Ny,kgrid.Nz);
for i = 1:kgrid.Nx
    for j = 1:kgrid.Ny
        if ~isnan(hologram_surf(i,j))
            val_hologram = round(hologram_surf(i,j));
            while val_hologram>0
            hologram_map(i,j,val_hologram)=1;
            val_hologram = val_hologram - 1;
            end
            val_t = round(mean(find(transducer(i,j,:)==1))); % remove the lens material behind the transducer 
            while val_t>0
            hologram_map(i,j,val_t)=0;
            val_t = val_t - 1;
            end
        end
    end
end

% set the medium parameters for the holographic lens 
medium.sound_speed(hologram_map~=0) = c_lens;
medium.density(hologram_map~=0)     = d_lens;
medium.alpha_coeff(hologram_map~=0) = alpha_lens;

% create the time array 
kgrid.t_array = 0:kgrid.dt:(kgrid.Nz*kgrid.dx*2)/c_water;
kgrid.Nt = length(kgrid.t_array);

% define the source as the transducer
source.p_mask = transducer;

% record the minimum pressure at every grid point
sensor.mask = ones(kgrid.Nx,kgrid.Ny,kgrid.Nz);
sensor.record = {'p_min'}; 

save('hologram_map.mat','hologram_map','-v7.3')
save('transducer.mat','transducer','-v7.3')

% run the simulation
input_args = {'PlotPML', true,'PMLInside', false, 'SaveToDisk', 'kwave_input_data_testgreytriple.h5'};
kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

sensor_data = h5read('kwave_output_data_testgreytriple.h5', '/p_min');

%% Visualisation
% Load test simulation data and normalize it
Data_norm = normalize((abs(sensor_data)),"range");
Data = reshape(Data_norm,kgrid.Nx,kgrid.Ny,kgrid.Nz);

% XY visualisation
Image_plane_data_XY = Data(:,:,focaldepth);
Image_plane_data_XY = normalize(reshape(Image_plane_data_XY,1,[]),"range");
Image_plane_data_XY = reshape(Image_plane_data_XY,kgrid.Nx,kgrid.Ny);

% display the simulated data
figure
imagesc(Image_plane_data_XY)
xticklabels = round(linspace(-round(kgrid.x_size*10^3/2),round(kgrid.x_size*10^3/2),5));
xticks = linspace(1, size(Image_plane_data_XY, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = round(linspace(round(kgrid.y_size*10^3/2),-round(kgrid.y_size*10^3/2),6));
yticks = linspace(1, size(Image_plane_data_XY, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
colormap(hot)
ylabel('y (mm)')
xlabel('x (mm)')
title('PU Simulation pressure map in the X-Y plane')
a=colorbar;
a.Label.String = 'Normalised pressure';
daspect([1 1 1])
clim([0 1])

% XZ visualisation 
Image_plane_data_XZ = squeeze(Data(round(kgrid.Nx/2),:,:));

% highlight the transducer and hologram
Image_plane_data_XZ(transducer(round(kgrid.Nx/2),:,:)==1)=1;
Image_plane_data_XZ(hologram_map(round(kgrid.Nx/2),:,:)==1)=1;

% display the target data
figure
imagesc(Image_plane_data_XZ)
xticklabels = round(linspace(0,kgrid.z_size*10^3,5));
xticks = linspace(1, size(Image_plane_data_XZ, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = round(linspace((kgrid.x_size/2)*10^3,-(kgrid.x_size/2)*10^3,5));
yticks = linspace(1, size(Image_plane_data_XZ, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
colormap(hot)
ylabel('x (mm)')
xlabel('z (mm)')
title('Flexible Resin Simulated pressure map in the X-Z plane')
a=colorbar;
a.Label.String = 'Normalised pressure';
daspect([1 1 1])
clim([0 1])

% YZ visualisation
Image_plane_data_YZ = squeeze(Data(:,round(kgrid.Ny/2),:));

% highlight the transducer and hologram
Image_plane_data_YZ(transducer(:,round(kgrid.Ny/2),:)==1)=1;
Image_plane_data_YZ(hologram_map(:,round(kgrid.Ny/2),:)==1)=1;

% display the target data
figure
imagesc(Image_plane_data_YZ)
xticklabels = round(linspace(0,kgrid.z_size*10^3,5));
xticks = linspace(1, size(Image_plane_data_YZ, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = round(linspace(-(kgrid.y_size/2)*10^3,(kgrid.y_size/2)*10^3,5));
yticks = linspace(1, size(Image_plane_data_XY, 2), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
colormap(hot)
ylabel('y (mm)')
xlabel('z (mm)')
title('Flexible Resin Simulated pressure map in the Y-Z plane')
a=colorbar;
a.Label.String = 'Normalised pressure';
daspect([1 1 1])
clim([0 1])
