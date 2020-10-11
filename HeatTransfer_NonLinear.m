%% 1D Nonlinear FEM Code for Heat Equation
%
% -div(k grad U) = f
%
% where the material model is non-linear i.e.
% k = k0(1+Beta U)
%
% This example is from Reddy Nonlinear FEM book.
%
% Author: Abdullah Waseem       
% Created: 04-November-2018
% Contact: engineerabdullah@ymail.com


% Initializing
clear; clc; clf; path(pathdef);
addpath FECore/

%% Creating 1D mesh
xstart = 0;     % start point
xend   = 1;     % End point
tne = 20;       % Total number of elements.

% Element type:         Q1 --> LINEAR,  Q2 --> QUADRATIC
elementtype = 'Q2';
% Creating 1D Mesh.
[ L, lnn, nne, el, egnn, tnn, x ] = CreateMesh( elementtype, tne, xstart, xend  );

%% Material Properties
k0   = 2e-1 * ones(tne,1);          % Conductivity
Beta = 2e-3 * ones(tne,1);          % Nonlinearity in the model

%% Finite Element Data
% Gauss Quadrature
ngp = 3;
run('GaussianLegendre');
% Shape Functions
run('ShapeFunctions')

%%	Boundary Conditions
tn = 1 : tnn;               % Node iterator
p = [tn(1) tn(end)];        % Precribed nodes
f = setdiff(tn, p);         % Free nodes

% Initializing solution, and solution at current Newton iteration
un = 273*ones(tnn,1);       

% Drichlet BC
un(1)   = 500;      
un(end) = 300;

up = un;                % solution at previous Newton iteration

%% Newton Raphson Algorithm

% Parameters
tol = 1e-8;             % Tolerance in the solution
maxiter = 20;           % Maximum number of iterations

% Initial calculation of Tangent Matrix 'T' and
% Internal Flux Vector 'FINT'
run('TangentFint')

% Neumann BC 
FEXT = zeros(tnn,1);

% Initializing 'du' vector
du   = zeros(tnn,1);

cnv = 1;        % Initializing convergence variable
iter = 1;       % Initializing iteration variable

% Newton Raphson Routine for Nonlinear Solution
while cnv > tol && iter<maxiter
    
	% Solving
	du(f,1) = T(f,f) \ ( FEXT(f,1) - FINT(f,1) );
    
    % Incrementing the solution vector
	un  = un + du;
	
    % Tangent Matrix 'T' and the Internal Flux Vector 'FINT'
	run('TangentFint')
	
    % Checking Convergence
	cnv = norm(up-un);
	
	disp(['               Newton Itr.: ' num2str(iter) '       Conv.:' num2str(cnv)]);
	
    % Assigning solution 'un', at curernt Newton iteration to the previous 'up' one
	up = un;
    
    % Incrementing iteration variable.
	iter=iter+1;
    
    % Plotting the solution
	plot(x,un,'LineWidth',2)
	drawnow
    
end
