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
		Jcbn = B(gs,:)*x(egnn(en,:));
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
        Ce(:,:,en) = Ce(:,:,en) + N(gs,:)'      * c  * N(gs,:)      * glw(gs) * Jcbn;
        % Elemental Conductivity Matrix
        Ke(:,:,en) = Ke(:,:,en) + B(gs,:)'/Jcbn * lm      * B(gs,:)/Jcbn * glw(gs) * Jcbn;
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
Beta  = 1;

% INITIALIZING
U = zeros(2*tnn,1);         % Total Vector -- Displacement + Temperature
v = zeros(tnn,1);           % Velocities
u = zeros(tnn,1);           % Displacements
w = 0*ones(tnn,1);          % Temperatures
U(tnn+1:2*tnn,1) = 0;
U(tnn+1,1) = 927; U(2*tnn,1) = 0;
F = zeros(2*tnn,1);

% ASSEMBLY OF THE BIG SYSTEM
A = [1/dt*barM+Beta^2*dt*barE -Beta*barG;
        Beta*barG'       1/dt*barC+Beta*barK];
Aff = decomposition(A(f,f));

% TIME STEPPING
En = zeros(tnts,1);
for t = 1 : tnts
    
    %Previous data
    up = u;
    vp = v;
    wp = U(tnn+1:2*tnn,1);
    
    % The total Forcing Vector
    F = [1/dt*barM*vp-barE*(up+Beta*(1-Beta)*dt*vp)+(1-Beta)*barG*wp;
         1/dt*barC*wp-(1-Beta)*barG'*vp-(1-Beta)*barK*wp];
    
    % Solving the Velocities and Temperature
    U(f,1) = Aff \ (F(f,1) - A(f,p)*U(p,1));
    
    %
    v = U(1:tnn,1);
    u = up + dt*(1-Beta)*vp + Beta*dt*v;
    w = U(tnn+1:2*tnn,1);
     
    figure(1); subplot(2,1,1,'align'); plot(x,u); title('Displacement');   drawnow
               subplot(2,1,2,'align'); plot(x,w); title('Temperature');    drawnow
           
    En(t,1) = 1/2*u'*barE*u + 1/2*v'*barM*v;
end

figure(1); subplot(2,1,1,'align'); plot(x,u); title('Displacement');   drawnow;
           subplot(2,1,2,'align'); plot(x,w); title('Temperature');    drawnow;

figure(2);
plot([1:tnts]'*dt,En)
title('Energy')


