% Kadir Ercan Özdemir
% Efficient Lennard-Jones potential energy calculation
% This version improves the efficiency of the original code by:
% 1) Removing duplicate pair calculations
% 2) Avoiding expensive square root operations
% 3) Using vectorized coordinate generation
% 4) Fixing the random seed for debugging

n = 30;     % Number of particles in the simulation
L = 10;     % Length of the cubic simulation box

% IMPORTANT CHANGE 1:
% Fix the random number generator seed so the same coordinates are
% generated every time the program runs. This is essential for debugging
% because it allows direct comparison between different versions of the code.
rng(1);

% IMPORTANT CHANGE 2:
% Instead of using a loop to assign coordinates one particle at a time,
% MATLAB can generate all coordinates at once using vectorized operations.
% Vectorized operations are much faster in MATLAB because they use
% optimized internal numerical routines.

xcoord = rand(n,1)*L;
ycoord = rand(n,1)*L;
zcoord = rand(n,1)*L;

energy = 0;   % Initialize total Lennard-Jones energy

% IMPORTANT CHANGE 3 (MAJOR EFFICIENCY IMPROVEMENT):
% The original code calculated each particle pair interaction twice:
% once when particle i interacted with particle j and again when
% particle j interacted with particle i.
%
% In the revised code, interactions are calculated only for j = i+1:n.
% This guarantees that every particle pair is evaluated exactly once,
% which reduces the number of calculations by roughly a factor of two.

for i = 1:n
    for j = i+1:n
        
        % Compute distance components between particle i and particle j
        rx = xcoord(j) - xcoord(i);
        ry = ycoord(j) - ycoord(i);
        rz = zcoord(j) - zcoord(i);
        
        % IMPORTANT CHANGE 4:
        % Instead of computing the distance using a square root
        % r = sqrt(rx^2 + ry^2 + rz^2)
        %
        % we compute the squared distance r^2. The Lennard-Jones
        % potential only requires inverse powers of r, so the
        % square root operation can be avoided. Square root
        % calculations are computationally expensive, so avoiding
        % them improves efficiency.

        r2 = rx^2 + ry^2 + rz^2;   % squared distance
        
        % Compute inverse powers of the distance using r^2
        inv_r2 = 1 / r2;
        inv_r6 = inv_r2^3;
        inv_r12 = inv_r6^2;
        
        % Lennard-Jones potential energy
        deltaE = 4*(inv_r12 - inv_r6);
        
        % IMPORTANT CHANGE 5:
        % In the original code, deltaE was divided by 2 because each
        % interaction was counted twice. Since each pair is now
        % evaluated only once, the division by 2 is no longer required.

        energy = energy + deltaE;
        
    end
end

% Display total Lennard-Jones potential energy
energy