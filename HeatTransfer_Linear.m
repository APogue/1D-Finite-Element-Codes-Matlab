%% 1D Nonlinear FEM Code for Transient Heat Equation
% 
% -div(k0 grad U) = 0
%
% Author: Abdullah Waseem       
% Created: 27-August-2017
% Contact: engineerabdullah@ymail.com

clear; clc; clf; path(pathdef);
addpath FECore/

%% 1D Meshing
xstart = 0;             % Start point
xend   = 1;             % End point
tne = 10;              % Total number of element in the domain.

% Element type:         Q1 --> LINEAR,  Q2 --> QUADRATIC
elementtype = 'Q1';
% Creating 1D Mesh.
[ L, lnn, nne, el, egnn, tnn, x ] = CreateMesh( elementtype, tne, xstart, xend  );
        
%% Material Properties (Constant with in elements -- Q0)
k0 = 100 * ones(tne,1);     % Conductivity
rho    = 1 * ones(tne,1);       % Density
c      = 1 * ones(tne,1);       % Heat capacity

% % Volume fraction of inclusion
% vf = 0.4;
% % Inlcusion material properties.
%   k0((tne/2) - vf*tne/2 : (tne/2) + vf*tne/2, 1) = 100*k0(1);
%  rho((tne/2) - vf*tne/2 : (tne/2) + vf*tne/2, 1) = 100*rho(1);
%    c((tne/2) - vf*tne/2 : (tne/2) + vf*tne/2, 1) = 100*c(1);

%% Finite Element Data
% Gauss Quadrature
ngp = 3;
run('GaussianLegendre.m');
% Shape Functions
run('ShapeFunctions.m');

%% 1D FEM CORE
% Initializing elemental Conductivity, Heat Capacity and Source Vector
Ke = zeros(nne, nne, tne);        
Ce = zeros(nne, nne, tne);        
Fe = zeros(nne,1,tne);
% Element loop
for en = 1 : tne
    % Gauss integration loop
	for gs = 1 : ngp
		
		% Jacobian Matrix
		Jcbn = B(gs,:)*x(egnn(en,:));
		
		% Iso-parameteric map
		x_z  = N(gs,:) * x(egnn(en,:));
		
		%Source at that gauss point
		source = x_z*sin(2*pi*x_z);  % This is an example
		
		% Element Conductivity Matrix
		Ke(:,:,en) = Ke(:,:,en) + B(gs,:)'/Jcbn * k0(en) * B(gs,:)/Jcbn * glw(gs) * Jcbn;
		
		% Element Heat Capacity Matrix
		Ce(:,:,en) = Ce(:,:,en) + N(gs,:)' * rho(en)*c(en) * N(gs,:) * glw(gs) * Jcbn;
		
		% Element Source Vector
		Fe(:,1,en) = Fe(:,1,en) + N(gs,:)' * source * glw(gs) * Jcbn;
	end
end
% Assemble barK, barC and barF
[ barK, barC, barF ] = Assembler( egnn, nne, tne, tnn, Ke, Ce, Fe, 'sparse' );

%% Apply Drichlet Boundary Condition
% FOR 1D STEADY STATE EXAMPLE

tn = 1 : tnn;               % Node iterator
p = [tn(1) tn(end)];                % Precribed nodes
f = setdiff(tn, p);         % Free nodes

% Initializing the solution vector. 
u0 = 273;
u = u0*ones(tnn,1); 

% Prescribed value of solution
u(1) = u0 + 0;
u(end) = u0 + 0;

% Solution.
u(f) = barK(f,f) \ ( - barK(f,p) * u(p) + barF(f,1) );

% Reaction flux
Fext = barK(p,p)*u(p) + barK(p,f)*u(f);

%% Plotting the results.
plot(x, u, 'Color', 'k', 'LineWidth', 1.2);




