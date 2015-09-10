% This script takes in the data for a PSD, and randomly generates particles
% within a designated area to approximate that PSD.
clear;
clc;
%tic;

%========================================
% USER INPUT DATA (IN SI UNITS)
%========================================

% The density of the simulated spheres
ball_density = 2650; % kg/m3

% x/y/z is the length/width/height of the box which no longer needs to be a
% perfect cube. As a general guideline, each box dimension should be >
% 10*{diameter of the largest particle}.
x = 1.0/1000; % m
y = 1.0/1000; % m
z = 1.0/1000; % m

% A text file is required containing the PSD by NUMBER (Q0) for the grading
% of interest. This file has 3 columns: lower class limit, upper class
% limit, percentage in class.
psd_filename = 'Grading_ToyouraSand_q0.txt'; % Currently limits are in mm in this file

% Optionally modify the PSD to disregard particles below some minimum size.
min_size = 0.00012; % m

% The void ratio for the sample. This controls the number of particles
% which will be placed inside the cell. The smaller the void ratio, the
% longer the script requires to run. Using a void ratio of around 1 is
% usually fairly quick.
void_ratio = 0.666; % Volume of voids / volume of solids

%========================================

% Take in the PSD by number (Q0) from a text file.
actual_psd = load(psd_filename);
actual_psd(:,1:2) = actual_psd(:,1:2)./1000; % Convert from mm to m
mean_size = 0.5.*(actual_psd(:,1) + actual_psd(:,2));

% Modify the PSD to disregard particles below some minimum size. Although
% the numbers of particles are not prohibitively high, the problem would
% arise when calculating the timestep, as it is directly proportional to
% the root of mass, so the presence of very small particles would
% necessitate using a tiny timestep.
for i = 1:length(mean_size)
    if actual_psd(i,2) <= min_size
        actual_psd(i,3) = 0;
    elseif (actual_psd(i,1) < min_size) && (actual_psd(i,2) > min_size)
        actual_psd(i,3) = actual_psd(i,3)*(actual_psd(i,2)-min_size)/(actual_psd(i,2) - actual_psd(i,1));
        actual_psd(i,1) = min_size;
    end
end

% Rescale the percentage categories
summation = sum(actual_psd(:,3));
actual_psd(:,3) = 100*actual_psd(:,3)./summation;

%========================================
% Now determine how many particles should be generated.
particle_volume = (x*y*z)/(1+void_ratio);
volume_of_100_particles = (1/6)*pi*sum((mean_size(:).^3).*actual_psd(:,3));
approx_no_particles_reqd = round(100*particle_volume/volume_of_100_particles);

fprintf(1,'The number of particles to be generated is approximately %i.\n',approx_no_particles_reqd);
pause(4); % Wait a few seconds

%========================================
% Draw plots of Q0 and Q3, to see what they look like
% vol_dist = 100.*(1/6)*pi*(mean_size(:).^3).*actual_psd(:,3)./volume_of_100_particles;
% 
% figure;
% hold on;
% plot(1000*mean_size,cumsum(actual_psd(:,3)),'b');
% plot(1000*mean_size,cumsum(vol_dist),'r');
% xlabel('Particle Size (mm)');
% ylabel('Size Distribution');
% legend('Number Size Distribution, Q0','Volume/Mass Size Distribution, Q3');
% legend boxoff;
% hold off;

%========================================
% Initially generate all of the particles in a list before placing them
% inside a box. It is most accurate to generate the particles by class, and
% ensure that the volumes in each class correspond as closely as possible
% to the required volume in that class.

% class_volume contains the particle volumes subdivided by class.
class_volumes = (1/6)*pi*(mean_size(:).^3).*actual_psd(:,3);
class_volume = particle_volume*class_volumes./sum(class_volumes);

% data is used to contain the diameters and coordinates of the generated
% particles. Initialise it at a size which is far larger than necessary to
% avoid having to extend the array later.
data = zeros(2*approx_no_particles_reqd,4);

% For each size class, starting at the largest, successively generate
% particles until the cumulative volume within the class is greater than
% that required (in class_volume). Keep all n particles if their total volume
% is closer to the required volume than the volume of (n-1) particles, or
% else delete the nth particle.

% Then calculate the volume disparity (positive if the total volume in the
% class is < the required volume, else negative). This volume is added to
% the volume in the next (smaller) size class, i.e., the volume disparities
% are "trickled down" to successively smaller size classes. Then move on to
% the next size class using the revised volume.

j = 1; % The row number in data

for i = length(actual_psd):-1:1
    accumulated_volume = 0; % Initialise at 0
    
    if class_volume(i) > 0
        while accumulated_volume < class_volume(i);
            data(j,1) = actual_psd(i,1) + rand*(actual_psd(i,2) - actual_psd(i,1));
            added_volume = (1/6)*pi*(data(j,1))^3;
            accumulated_volume = accumulated_volume + added_volume;
            j = j + 1;
        end
        
        if (accumulated_volume - class_volume(i)) > (class_volume(i) - (accumulated_volume - added_volume))
            % Delete element j - it will actually be overwritten by the particles
            % in the next class rather than deleted, except for the final
            % class which will be deleted.
            j = j - 1;
        end
    end
end

% Note that the final class (i == 1) must be treated as a special case to
% decrement the j counter.
j = j - 1;

%========================================
% Sort the data array by descending particle diameter and delete the
% redundant rows of zeros at the end of the array
data = sort(data,'descend');

data(j+1:end,:) = [];

vol_spec = (1/6)*pi*sum((data(:,1).^3));
fprintf(1,'\nExactly %i particles were generated: their volume is within %1.2f%% of the required volume.\n\n',j,100*abs(particle_volume-vol_spec)/particle_volume);
pause(4); % Wait a few seconds

%========================================
% Create the data for plotting a PSD
%diams = flipud(data(:,1)); % Sort from small to large
%cumvolumes = 100*cumsum((pi/6)*(diams.^3))/vol_spec;
%diams = 1000*diams; % Convert to mm before saving
%save('1.0_Generated_PSD.mat', 'diams', 'cumvolumes'); % Save these to .mat file

%========================================
% Create a tiny tolerance equivalent to 0.00000001*{minimum particle
% diameter}.
tol = 0.00000001*data(j,1);

% Now sequentially place the particles inside a confining volume, starting
% wih the largest particle and finishing with the smallest. The particles
% are not permitted to overlap. It is most efficient to subdivide the large
% volume into many smaller subvolumes, say 6^3, and maintain separate lists
% for each subvolume.

% Divide into 216 (6^3) smaller cuboids.
lcube = 1;
ncubes = lcube^3;

% The structure will contain the diameter of the particles and the
% coordinates of placed particles.
data2 = struct([]);

% It is important to note the way in which data2 is organised:
% data.(cell_id).diameter gives the diameters of all particles in the list
% cell_id. Similarly for x, y and z separately.

cube_bounds = zeros(ncubes,6);

% The following are defined to save on calculations within the for loop
% below.
radii = 0.5*data(:,1);
ncoverlc = ncubes/lcube;

% Set the dimensions for each cube.
for i = 1:ncubes
    % First change the x direction, then y and finally z
    cube_bounds(i,1) = mod(i-1,lcube)*x/lcube - (radii(1) + tol);
    cube_bounds(i,3) = mod(floor((i-1)/lcube),lcube)*y/lcube - (radii(1) + tol);
    cube_bounds(i,5) = (floor((i-1)/ncoverlc))*z/lcube - (radii(1) + tol);
    
    % Create the list of cell names
    allnames{i} = strcat('id_',num2str(i));
end

cube_bounds(:,2) = cube_bounds(:,1) + x/lcube + data(1,1) + tol;
cube_bounds(:,4) = cube_bounds(:,3) + y/lcube + data(1,1) + tol;
cube_bounds(:,6) = cube_bounds(:,5) + z/lcube + data(1,1) + tol;

maxtries = 500000000; % The maximum number of attempts to make at placing an atom
celldims = [x; y; z]; % The dimensions of the periodic cell

% This loop cycles through each placed particle to ensure no overlaps. The
% first particle is placed trivially.
for n = 1:j
    placed = 0;
    tries = 0; % The number of attempts made at placing an atom
    
    tic
    while (placed == 0) && (tries < maxtries)
        coord = celldims.*rand(3,1); % Randomly select a particle position
        
        placed = 1;
        tries = tries + 1;
        memberlist = []; % Empty this list of cube ids
        sc = []; % And this which is for reducing the number of entries in the for loop below
        
        % Find the cube(s) this belongs to. Note that the z direction has
        % the greatest runs of unchanged bounds. Therefore, it would be
        % most efficient to begin with that.
        for i = 1:lcube
            if (coord(3) > cube_bounds(1+(i-1)*ncoverlc,5)) && (coord(3) < cube_bounds(1+(i-1)*ncoverlc,6))
                sc((length(sc)+1):(length(sc)+ncoverlc)) = (1+(i-1)*ncoverlc):(i*ncoverlc);
            end
        end
        
        for i = 1:length(sc)
            if (coord(1) > cube_bounds(sc(i),1) && coord(1) < cube_bounds(sc(i),2)) && (coord(2) > cube_bounds(sc(i),3) && coord(2) < cube_bounds(sc(i),4))
                % If it is a member
                
                name = allnames{sc(i)};
                
                if isfield(data2,name) == 1
                    % These are overwritten each time, so do not require
                    % re-initialisation or deletion
                    dist = (((coord(1) - data2.(name).x).^2 + (coord(2) - data2.(name).y).^2 + (coord(3) - data2.(name).z).^2).^0.5)';
                    gap = dist - (radii(data2.(name).index) + radii(n) + tol);
                    
                    if min(gap) < 0
                        placed = 0;
                        break;
                    end
                end
                
                memberlist(length(memberlist)+1) = sc(i);
            end
        end
        
        if placed == 1
            % If the particle protrudes into the "squish" volumes at the edge of
            % the periodic boundaries, create a copy of the particle with
            % mapped coordinates and add this to the appropriate membership lists
            periodic = zeros(6,1);
            for i = 1:3
                if (coord(i)-radii(n)) < 0
                    periodic(2*i-1) = 1;
                elseif (coord(i)+radii(n)) > celldims(i)
                    periodic(2*i) = 1;
                end
            end
            
            if sum(periodic) > 0 % If there is a protruding particle
                second_coord = {}; % Initialise here
                second_coord{1} = coord;
                counter = 2;
                
                for i = 1:6
                    if (periodic(i) == 1) && (mod(i,2) == 1)
                        if i < 2
                            second_coord{counter} = [coord(1)+x coord(2) coord(3)];
                        elseif i < 4
                            second_coord{counter} = [coord(1) coord(2)+y coord(3)];
                        else
                            second_coord{counter} = [coord(1) coord(2) coord(3)+z];
                        end
                        counter = counter + 1;
                    elseif (periodic(i) == 1) && (mod(i,2) == 0)
                        if i == 2
                            second_coord{counter} = [coord(1)-x coord(2) coord(3)];
                        elseif i == 4
                            second_coord{counter} = [coord(1) coord(2)-y coord(3)];
                        else
                            second_coord{counter} = [coord(1) coord(2) coord(3)-z];
                        end
                        counter = counter + 1;
                    end
                end
                
                % A particle could potentially need to be mapped to more than three
                % positions; a particle which has a centrepoint on the corner of
                % the assembly would need to be mapped to all corners.
                if sum(periodic) > 1
                    minimum = 10*celldims; % Initialise at large values
                    maximum = -10*celldims;
                    
                    for i = 1:length(second_coord)
                        for q = 1:3
                            minimum(q) = min(second_coord{i}(q),minimum(q));
                            maximum(q) = max(second_coord{i}(q),maximum(q));
                        end
                    end
                    
                    diff = min(abs(maximum - minimum));
                    
                    if diff > 0 % This means that none of the entries in maximum are the same as the equivalent in minimum
                        % Need all 8 mappings; already have 4 (including original)
                        second_coord{counter} = [coord(1) minimum(2)+maximum(2)-coord(2) minimum(3)+maximum(3)-coord(3)];
                        second_coord{counter+1} = [minimum(1)+maximum(1)-coord(1) coord(2) minimum(3)+maximum(3)-coord(3)];
                        second_coord{counter+2} = [minimum(1)+maximum(1)-coord(1) minimum(2)+maximum(2)-coord(2) coord(3)];
                        second_coord{counter+3} = [minimum(1)+maximum(1)-coord(1) minimum(2)+maximum(2)-coord(2) minimum(3)+maximum(3)-coord(3)];
                    else
                        % Need 4 mappings, have 3
                        second_coord{counter} = [minimum(1)+maximum(1)-coord(1) minimum(2)+maximum(2)-coord(2) minimum(3)+maximum(3)-coord(3)];
                    end
                end
                
                % Now it is necessary to determine if overlaps exist for these sets
                % of coordinates. This is done in the same way as before. Exclude
                % the first (real) coordinate as this was dealt with
                % separately.
                for q = 2:length(second_coord)
                    sc = [];
                    
                    for i = 1:lcube
                        if (second_coord{q}(3) > cube_bounds(1+(i-1)*ncoverlc,5)) && (second_coord{q}(3) < cube_bounds(1+(i-1)*ncoverlc,6))
                            sc((length(sc)+1):(length(sc)+ncoverlc)) = (1+(i-1)*ncoverlc):(i*ncoverlc);
                        end
                    end
                    
                    for i = 1:length(sc)
                        if (second_coord{q}(1) > cube_bounds(sc(i),1) && second_coord{q}(1) < cube_bounds(sc(i),2)) && (second_coord{q}(2) > cube_bounds(sc(i),3) && second_coord{q}(2) < cube_bounds(sc(i),4))
                            name = ['id_',num2str(sc(i))];
                            
                            if isfield(data2,name) == 1
                                % These are overwritten each time, so do not require
                                % re-initialisation or deletion
                                dist = (((second_coord{q}(1) - data2.(name).x).^2 + (second_coord{q}(2) - data2.(name).y).^2 + (second_coord{q}(3) - data2.(name).z).^2).^0.5)';
                                gap = dist - (radii(data2.(name).index) + radii(n) + tol);
                                
                                if min(gap) < 0
                                    placed = 0;
                                    break;
                                end
                            end
                        end
                    end
                    
                    if placed == 0
                        break;
                    end
                end
            end
        end
    end
    
    if (tries == maxtries) && (placed == 0)
        break;
    end
    
    % Write the coordinates of the newly-placed particle to the global file
    % and also to the cube sub-lists
    data(n,2:4) = coord;
    
    for i = 1:length(memberlist)
        name = ['id_',num2str(memberlist(i))];
        
        if isfield(data2,name) == 1
            leng = length(data2(1).(name).index);
        else
            leng = 0;
        end
        
        data2(1).(name).index(leng+1) = n;
        data2(1).(name).x(leng+1) = coord(1);
        data2(1).(name).y(leng+1) = coord(2);
        data2(1).(name).z(leng+1) = coord(3);
    end
    
    if sum(periodic) > 0
        for q = 2:length(second_coord)
            sc = [];
            
            for i = 1:lcube
                if (second_coord{q}(3) > cube_bounds(1+(i-1)*ncoverlc,5)) && (second_coord{q}(3) < cube_bounds(1+(i-1)*ncoverlc,6))
                    sc((length(sc)+1):(length(sc)+ncoverlc)) = (1+(i-1)*ncoverlc):(i*ncoverlc);
                end
            end
            
            for i = 1:length(sc)
                if (second_coord{q}(1) > cube_bounds(sc(i),1) && second_coord{q}(1) < cube_bounds(sc(i),2)) && (second_coord{q}(2) > cube_bounds(sc(i),3) && second_coord{q}(2) < cube_bounds(sc(i),4))
                    name = ['id_',num2str(sc(i))];
                    
                    if isfield(data2,name) == 1
                        leng = length(data2(1).(name).index);
                    else
                        leng = 0;
                    end
                    
                    data2(1).(name).index(leng+1) = n;
                    data2(1).(name).x(leng+1) = second_coord{q}(1);
                    data2(1).(name).y(leng+1) = second_coord{q}(2);
                    data2(1).(name).z(leng+1) = second_coord{q}(3);
                end
            end
        end
    end
    
    toc
    fprintf(1,'Placed particle %i of %i in %i attempts.\n',n,j,tries); % To keep an eye on progress
    fprintf(1, 'wall time: %f %f',toc, toc/tries);
end

if (tries == maxtries) && (placed == 0)
    fprintf(1,'\nPlaced only %i of %i required particles.\n',n,j);
else
    fprintf(1,'\n\nPlaced all %i required particles.\n',j);
end

%========================================
% Write the particle coordinates to a text file, openable using the read_data command in LAMMPS.
filename = strcat('Particle_Coordinates_',num2str(void_ratio),'.lj');
fprintf(1,'Writing data to the text file ''%s''....\n',filename);
fid = fopen(filename,'wt');
fprintf(fid,'%s\n',['# Particle coordinates for ' num2str(j) ' spheres with diameters between ' num2str(data(j,1)) ' and ' num2str(data(1,1)) ' m.']);
fprintf(fid,'%s\n','1 atom types');
fprintf(fid,'%s\n',strcat(num2str(j),' atoms'));

% Optionally shift the positions of the particles inside the periodic cell.
% If all of these numbers are zero, the particles will not be shifted.
xlo = 0.0;
ylo = 0.0;
zlo = 0.0;

xhi = xlo + x;
yhi = ylo + y;
zhi = zlo + z;

fprintf(fid,'%1.16e %1.16e %s\n',xlo,xhi,'xlo xhi');
fprintf(fid,'%1.16e %1.16e %s\n',ylo,yhi,'ylo yhi');
fprintf(fid,'%1.16e %1.16e %s\n',zlo,zhi,'zlo zhi');
fprintf(fid,'%s\n',' ');
fprintf(fid,'%s\n','Atoms');
fprintf(fid,'%s\n',' ');
for m = 1:j
    fprintf(fid,'%i %i %1.16e %1.16e %1.16e %1.16e %1.16e\n',m,1,data(m,1),ball_density,data(m,2)+xlo,data(m,3)+ylo,data(m,4)+zlo);
end
fprintf(fid,'%s\n',' ');
fprintf(fid,'%s\n','Velocities');
fprintf(fid,'%s\n',' ');
% All particles are given an initial velocity of zero
for m = 1:j
    fprintf(fid,'%i %f %f %f %f %f %f\n',m,0.0,0.0,0.0,0.0,0.0,0.0);
end
fclose(fid);
fprintf(1,'Finished writing data!\n\n');

%toc;
return
