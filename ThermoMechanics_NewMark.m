%% 1-D COUPLED THERMO-MECHANICS FINITE ELEMENT CODE
% 
% Equations: \nabla\cdot\sigma + \rho a = 0
% 	     \nabla\cdot j     + \rho c \dot{\theta} = 0
% 	     \sigma: stess, rho: density, j: heat flux, c: heat capacity
% 	     a: acceleration, \theta: temperature
%
% Author: Abdullah Waseem       
% Created: 19-March-2019
% Contact: engineerabdullah@ymail.com

clear; clc; clf; path(pathdef); format long
addpath FECore/

%% 1D Meshing
xstart = 0;             % Start point
xend   = 3;             % End point
tne    = 60;           % Total number of element in the domain.

% Element type:         Q1 --> LINEAR,  Q2 --> QUADRATIC
elementtype = 'Q2';
% Creating 1D Mesh.
[ L, lnn, nne, el, egnn, tnn, x ] = CreateMesh( elementtype, tne, xstart, xend  );
        
%% Material Properties (Constant with in elements -- Q0)

% MECHANICAL
E     = 2.1e8;      % Elasticity Tensor
Alpha = 11e-6;      % Coefficient of Thermal Expansion
rho   = 7850;      % Density

% THERMAL
Tref   = 273;     % Reference Temperature
c      = 500;     % Heat Capacity
lm     = .0043;     % Thermal Conductivity

%% Pre-calculation of Gauss-Legendre Quadrature, Shape function and their Derivatives
% Gauss Quadrature
ngp = 3;
run('GaussianLegendre.m');
% Shape Functions
run('ShapeFunctions.m');

%% 1D FEM CORE

% Initializing Element Matrices
Me = zeros(nne, nne, tne);      % Mass
Ee = zeros(nne, nne, tne);      % Stiffness
Ge = zeros(nne, nne, tne);      % Coupling
Ce = zeros(nne, nne, tne);      % Capacity
Ke = zeros(nne, nne, tne);      % Conductivity
Fe = zeros(nne, 1  , tne);      % Force

% Element loop
for en = 1 : tne
    % Gauss integration loop
	for gs = 1 : ngp
		
		% Jacobian Matrix
		Jcbn = B(gs,:) * x(egnn(en,:));
		% Iso-parameteric map
		x_z  = N(gs,:) * x(egnn(en,:));
		%Force at that gauss point
		force = (3*x_z + x_z^2)*exp(x_z);  % This is an example
		
		% Element Mass Matrix
		Me(:,:,en) = Me(:,:,en) + N(gs,:)'      * rho     * N(gs,:)      * glw(gs) * Jcbn;
        % Element Stiffness Matrix
		Ee(:,:,en) = Ee(:,:,en) + B(gs,:)'/Jcbn * E       * B(gs,:)/Jcbn * glw(gs) * Jcbn;
        % Element Coupling Matrix
        Ge(:,:,en) = Ge(:,:,en) + B(gs,:)'/Jcbn * Alpha*E * N(gs,:)      * glw(gs) * Jcbn;
        % Elemental Capacity Matrix
        Ce(:,:,en) = Ce(:,:,en) + N(gs,:)'      * c/Tref  * N(gs,:)      * glw(gs) * Jcbn;
        % Elemental Conductivity Matrix
        Ke(:,:,en) = Ke(:,:,en) + B(gs,:)'/Jcbn * lm/Tref * B(gs,:)/Jcbn * glw(gs) * Jcbn;
		% Element Force Vector
		Fe(:,1,en) = Fe(:,1,en) + N(gs,:)' * force * glw(gs) * Jcbn;
        
	end
end

% Assemble barK, barC and barF
[ barM, barE, barF ] = Assembler( egnn, nne, tne, tnn, Me, Ee, Fe, 'sparse' );
[ barG, barC, ~    ] = Assembler( egnn, nne, tne, tnn, Ge, Ce, Fe, 'sparse' );
[ barK, ~   , ~    ] = Assembler( egnn, nne, tne, tnn, Ke, Ke, Fe, 'sparse' );

%% BOUNDARY CONDITIONS

%  MECHANICAL/THERMAL -- FIXED AT BOTH ENDS         
p = [1 tnn 1+tnn 2*tnn];               % Prescribed 
f = setdiff(1:2*tnn,p);                % Free

%% NEWMARK + MONOLITHIC

% TIME DATA
T    = 600;                 % Total Time
dt   = 10;                 % Time Step Size
tnts = T/dt+1;              % Total Number of Time Steps

% TIME SCHEME
Beta  = 1/4;
Gamma = 1/2;

% INITIALIZING
u = zeros(tnn,1);           % Displacements
v = zeros(tnn,1);           % Velocities
a = zeros(tnn,1);           % Accelerations
w = zeros(tnn,1);           % Temperatures
w(1) = 927; w(end) = 0;
U = [u;w];                  % Total Vector -- Displacement + Temperature

% ASSEMBLY OF THE BIG SYSTEM
A = [barM/(Beta*dt^2)+barE -barG;
             barG' dt*barK+barC];
% When the system matrices does not change i.e. linear case
% The system matrix can be assembled, combined and decomposed for faster simulations
% This feature was first introduced in Matlab 2017b. If you have an older version of 
% Matlab then remove the word "decomposition" from the following line.
Aff = decomposition(A(f,f));

% TIME STEPPING
En = zeros(tnts,1);
for t = 1 : tnts
    %t
    %Previous data
    up = u;
    vp = v;
    ap = a;
    wp = w;
    
    % The total Forcing Vector
    F = [barM*up/(Beta*dt^2)-barM*vp/(Beta*dt)-(0.5-Beta)/Beta*ap;
         barC*wp+barG'*up];
    
    % Solving the Velocities and Temperature
    U(f,1) = Aff \ (F(f,1) - A(f,p)*U(p,1));
    
    %
    u = U(1:tnn,1);
    a = (u-up)/(Beta*dt^2) - vp/(Beta*dt) - (0.5-Beta)/Beta*ap;
    v = vp + (1-Gamma)*dt*ap + Gamma*dt*a;
    w = U(1+tnn:2*tnn,1);
     
    En(t) = 1/2*u'*barE*u + 1/2*v'*barM*v;
    
figure(1); subplot(2,1,1,'align'); plot(x,u); title('Displacement');   drawnow;
           subplot(2,1,2,'align'); plot(x,w); title('Temperature');    drawnow;
    
end

figure(1); subplot(2,1,1,'align'); plot(x,u); title('Displacement');   drawnow;
           subplot(2,1,2,'align'); plot(x,w); title('Temperature');    drawnow;


figure(2);
plot([1:tnts]'*dt,En);
title('Energy')


