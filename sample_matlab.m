% Sample matlab code for calculating LJ interactions without periodic
% boundary conditions, reduced units are used throughout
% Note that % indicates the beginning of a comment in the code and ;
% suppresses output to the screen 
% Set number of particles to include in calculation
n = 10;
% Define simulation box size (in reduced units)
L = 10;
% Define arrays for x, y, and z coords 
xcoord = zeros(n,1);
ycoord = zeros(n,1);
zcoord = zeros(n,1);
% Assign all coordinates random values within simulation box
for i = 1:n
    xcoord(i) = rand*L;
    ycoord(i) = rand*L;
    zcoord(i) = rand*L;
end
% Sum pairwise interactions between particles with LJ potential
energy = 0;
for i = 1:n
    for j = 1:1:i-1
        rx = xcoord(j)-xcoord(i);
        ry = ycoord(j)-ycoord(i);
        rz = zcoord(j)-zcoord(i);
        r = sqrt(rx*rx + ry*ry + rz*rz);
        deltaE = 4*(r^-12 - r^-6);
        energy = energy + deltaE/2;
    end
    for j = i+1:1:n
        rx = xcoord(j)-xcoord(i);
        ry = ycoord(j)-ycoord(i);
        rz = zcoord(j)-zcoord(i);
        r = sqrt(rx*rx + ry*ry + rz*rz);
        deltaE = 4*(r^-12 - r^-6);
        energy = energy + deltaE/2;
    end
end
%Print total energy to screen
energy

