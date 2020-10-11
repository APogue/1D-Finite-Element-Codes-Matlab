%% 1D Nonlinear FEM Code for Transient Heat Equation
% 
% -div(k0 grad U) + rho c \dot{U} = 0
%
% Author: Abdullah Waseem       
% Created: 27-August-2017
% Contact: engineerabdullah@ymail.com

clear; clc; close all; path(pathdef);
addpath FECore/

%% 1D Meshing
xstart = 0;             % Start point
xend   = 1;             % End point
tne = 100;              % Total number of element in the domain.

% Element type:         Q1 --> LINEAR,  Q2 --> QUADRATIC
elementtype = 'Q2';
% Creating 1D Mesh.
[ L, lnn, nne, el, egnn, tnn, x ] = CreateMesh( elementtype, tne, xstart, xend  );
        
%% Material Properties (Constant with in elements -- Q0)
k0 = 100 * ones(tne,1);     % Conductivity
rho    = 1 * ones(tne,1);       % Density
c      = 1 * ones(tne,1);       % Heat capacity

% Volume fraction of inclusion
vf = 0.4;
% Inlcusion material properties.
  k0((tne/2) - vf*tne/2 : (tne/2) + vf*tne/2, 1) = .1*k0(1);
 rho((tne/2) - vf*tne/2 : (tne/2) + vf*tne/2, 1) = 10*rho(1);
   c((tne/2) - vf*tne/2 : (tne/2) + vf*tne/2, 1) = 1*c(1);

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
		x_z  = 0;%N(gs,:) * x(egnn(en,:));
		
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

%% Transient Analysis.

nu = 1;                 % Parameter in generalized-theta time integration scheme
						% For Euler-forward use: nu = 0;
						% For Euler-backward use: nu = 1;
						% For Crank-Nicolson use: nu = 0.5 or 0.51

TotalTime = 1;          % Total simulation time
tnts = 100;             % Total number of time steps
dt = TotalTime/tnts;    % Time-step size
titrt = 0 : TotalTime/tnts : TotalTime; % Time iterator

%% Boundary Conditions
u0 = 273;               % Initial temperature
u_max = 1;              % Maximum tempearture

% Sinusoidal Loading.
osc = 1;                            % Number of oscillations
a = 2*pi*osc * titrt/TotalTime;     % Angle
u_p = u_max * sin(a) + u0;          % Temperature at the prescribed nodes.

% Initializing the solution
u  = u0 * ones(tnn,tnts+1);

tn = 1 : tnn;               % Node iterator
p = [tn(1)];                % Precribed nodes
f = setdiff(tn, p);         % Free nodes

% Partitioning the matrices
Kpp = barK(p,p); Kpf = barK(p,f); Kfp = barK(f,p); Kff = barK(f,f);
Cpp = barC(p,p); Cpf = barC(p,f); Cfp = barC(f,p); Cff = barC(f,f);

% When the conductivity and the capacity matrices does not change i.e. linear case
% The system matrix can be assembled, combined and decomposed for faster simulations
% This feature was first introduced in Matlab 2017b. If you have an older version of 
% Matlab then remove the word "decomposition" from the following line.
dA = decomposition( Cff + dt*nu*Kff );

% Initializing the reaction flux at the prescribed nodes.
Fext = zeros(2,tnts);

% Time integration loop
for t = 2 : tnts
	
	% Prescribing temperature in time.
	u(p,t)  = u_p(t);
	
    % Internal Flux Vector.
	Ft = dt*(1-nu)*Kfp*u(p,t-1) + dt*nu*Kfp*u(p,t) + dt*(1-nu)*Kff*u(f,t-1) ...
			+ Cfp*(u(p,t)-u(p,t-1)) - Cff*u(f,t-1);
    
    % Solution 
	u(f, t) = - dA  \  Ft;
    
    % Reaction flux at prescribed nodes
    Fext(:,t) = Kpp*u(p,t) + Kpf*u(f,t) + Cpp*(u(p,t) - u(p,t-1))/dt + Cpf*(u(f,t) - u(f,t-1))/dt;
    
    % Plotting the temperature field.
	plot(x,u(:,t)); ylim([u0-u_max u0+u_max]);
	title([num2str(t)])
	drawnow

end
