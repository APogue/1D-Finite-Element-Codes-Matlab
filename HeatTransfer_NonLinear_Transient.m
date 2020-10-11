%% 1D Nonlinear FEM Code for Transient Heat Equation
% 
% -div(k grad U) + rho \dot{U} = 0
%
% where the material model is non-linear i.e.
% k = k0(1+Beta U)
%
% This example is from Reddy Nonlinear FEM book.
%
% Author: Abdullah Waseem       
% Created: 04-November-2018
% Contact: engineerabdullah@ymail.com

clear; clc; clf; path(pathdef);
addpath FECore/

%% Creating 1D mesh
xstart = 0;         % start point
xend   = .5;        % End point
tne = 20;           % Total number of elements.

% Element type:         Q1 --> LINEAR,  Q2 --> QUADRATIC
elementtype = 'Q2';
% Creating 1D Mesh.
[ L, lnn, nne, el, egnn, tnn, x ] = CreateMesh( elementtype, tne, xstart, xend  );

%% Material Properties (Constant with in elements -- Q0)
k0   = 1e-1 * ones(tne,1);      % Conductivity
Beta = 1e-1 * ones(tne,1);      % Nonlinearity in the model
rho =  100 * ones(tne,1);       % Heat capacity per unit volume

% Volume fraction of inclusion
vf = 0.4;
% Inlcusion material properties.
  k0((tne/2) - vf*tne/2 : (tne/2) + vf*tne/2, 1) = 1e-1/2;
Beta((tne/2) - vf*tne/2 : (tne/2) + vf*tne/2, 1) = 1e-1/2;
 rho((tne/2) - vf*tne/2 : (tne/2) + vf*tne/2, 1) = 100/2;

%% Finite Element Data
% Gauss Quadrature
ngp = 3;
run('GaussianLegendre');
% Shape Functions
run('ShapeFunctions')

%% Assembly of solution and time independent FE matrix and vector.
% Initializing Mass and Force Vector
Me = zeros(nne, nne, tne);
Fe = zeros(nne,1,tne);

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
		Me(:,:,en) = Me(:,:,en) + N(gs,:)' * rho(en) * N(gs,:) * glw(gs) * Jcbn;
		
		% Element Force Vector
		Fe(:,1,en) = Fe(:,1,en) + N(gs,:)' * force * glw(gs) * Jcbn;
    end
end

%Assembly
[ M, ~, Fm ] = Assembler( egnn, nne, tne, tnn, Me, Me, Fe, 'sparse' );

%% Time Data.
TotalTime = 4;      
dt = 0.01;              % Time-step size
titrt = 0:dt:TotalTime; % Time iterator
tnts = size(titrt,2);   % Total number of time steps

%% Boundary Conditions
tn = 1 : tnn;               % Node iterator
p = [tn(1) tn(end)];        % Precribed nodes
f = setdiff(tn, p);         % Free nodes

% Initializing solution
iniu = 273;                 % Initial temperature
u = iniu*ones(tnn,tnts);

% Drichlet BC in time
u(1,:)   = 30*sin(2*pi*2*titrt/TotalTime)+273;
u(end,:) = 273;

%% Newton Raphson

% Parameters
tol = 1e-8;             % Tolerance in the solution
maxiter = 20;           % Maximum number of iterations

% Initial calculation of Tangent Matrix 'T' and
% Internal Force Vector 'FINT'
un = u(:,1);            % solution at current Newton iteration
up=un;                  % solution at previous Newton iteration
run('TangentFint')

% Neumann BC 
FEXT = zeros(tnn,1);

% Initializing 'du' vector
du   = zeros(tnn,1);

% Time integration loop
for t = 2 : tnts

	cnv = 1;        % Initializing convergence variable
	iter = 1;       % Initializing iteration variable
	
    % Assigning solution from previous time step to current
    un(p,1) = u(p,t);
	un(f,1) = up(f,1);
    
	disp(['At time: ' num2str(t*dt-dt) ' sec. and time step: ' num2str(t) ])
    
    % Newton Raphson Routine for Nonlinear Solution
	while cnv > tol && iter<maxiter
	
        % Solving
		du(f,1) = (M(f,f)+dt*T(f,f)) \ (-dt*FINT(f,1) - M(f,f)*(un(f,1)-u(f,t-1)) - M(f,p)*(un(p,1)-u(p,t-1)));
		
        % Incrementing the solution vector
		un = un + du;

        % Tangent Matrix 'T' and the Internal Flux Vector 'FINT'
		run('TangentFint')
        
        % Checking Convergence
		cnv = norm(up-un);
		
		disp(['               Newton Itr.: ' num2str(iter) '       Conv.:' num2str(cnv)]);
		
        % Assigning solution 'un', at curernt Newton iteration to the previous 'up' one
		up = un;
        
        % Incrementing iteration variable.
		iter=iter+1;
		
	end

	% Storing the converged solution in time. 
	u(:,t) = un;
    
    % Plotting the solution
	plot(x,u(:,t),'LineWidth',2); axis tight; ylim([240 306]); axis square;
	drawnow
	disp(' ')
	
end
