%% TRANSIENT DATA-DRIVEN SOLVER FOR DIFFUSION TYPE PROBLEMS (1D)
%  ------------
% Author: Abdullah Waseem       
% Created: 06-April-2020
% Contact: engineerabdullah@ymail.com
% 
% This finite element code solves a transient diffusion problem in a one-dimensional domain. The mass balance equation
% considered is " div(j) + dot(c) = 0 ". The chemical potential field is solved at the nodes while the concentrations are
% evaluated at the Gauss quadrature points. 

%% Initializing Code

clear; 
clc; 
close all;
format short;

path(pathdef);
addpath FECore/

%% 1D Meshing
xstart = 0;             % Start point
xend   = 10;            % End point
tne = 10;               % Total number of element in the domain.
% Element type:         Q1 --> LINEAR,  Q2 --> QUADRATIC
elementtype = 'Q1';
% Creating 1D Mesh.
[ L, lnn, nne, el, egnn, tnn, x ] = CreateMesh( elementtype, tne, xstart, xend  );

%% TEMPORAL DISCRETIZATION (for Backward-Euler)
T    = 100;             % Total Time
nTs  = 100;             % Number of time steps
dT   = T/nTs;           % Time Step Size
tItrt= 0:dT:T;          % Time iterator

%% Pre-calculation of Gauss-Legendre Quadrature, Shape function and their Derivatives
ngp = 2;                    % Gauss Quadrature
run('GaussianLegendre.m');

run('ShapeFunctions.m')     % Shape Functions

%% DATA-GENERATION

% REFERENCE MATERIAL PROPERTIES
Dref = 10;      % The reference diffusivity
nna=1;          % scaling diffusivity in the distance function
Lref = 1;       % The reference chemical modulus
nnb=1;          % scaling chemical modulus in the distance function
C0 = 1;         % Reference concentration
C1 = 3;         % Maximum concentration


nDP  = 20;              % Number of data points
DBSize = (nDP)*(nDP);   % Data-base size
mu1 = Lref*log(C1/C0);  % Maximum chemical potential value (it depends on what constitutive model you use)
GradRange = mu1/L;      % The gradient of the chemical potential (the data has to be consistent with the boundary conditions)

% The data-set is generated using the following constitutive model. 
% j = -Dref c gradMu
% Mu = Lref log (c/c0)
% where, j is the mass flux, c is the concentration, Mu is the chemical potential and gradMu is the gradient of the chemical
% potential.

% Concentrations
c_prm = zeros(1,DBSize);
for n = 1 : nDP
    c_prm((n-1)*nDP+1:n*nDP) = linspace(C0,C1,nDP);
end

% Chemical Potential
Mu_prm = Lref*log(c_prm/C0);

% Gradient of Chemical Potential
gradMu_prm = zeros(1,DBSize);
for n = 1 : nDP
    gradMu_prm((n-1)*nDP+1:n*nDP) = -GradRange+n*1*GradRange/(nDP);
end

% Mass Flux
tenJ_prm = -Dref*(c_prm).*gradMu_prm;

% ADDING NOISE TO THE GENERATED DATA-SET
    % to concentrations
    absErr = 1.0e-2*(C1-C0);
    relErr = 5.0e-3*(C1-C0);
    c_prm = c_prm + absErr*rand(1,DBSize) + relErr*rand(1,DBSize);
    % to chemical-potential
    absErr = 2.0e-2*Lref*log(C1/C0);
    relErr = 3.0e-2*Lref*log(C1/C0);
    Mu_prm = Mu_prm + absErr*rand(1,DBSize)+ relErr*rand(1,DBSize);
    % to mass flux
    absErr = 1.0e-2*C0*Dref(1,1)*GradRange;
    relErr = 3.0e-2*C0*Dref(1,1)*GradRange;
    tenJ_prm = tenJ_prm + absErr*rand(1,DBSize) + relErr*rand(1,DBSize);

% VISUALIZING THE DATA-SET
figure(1)
subplot(1,2,1)
plot3(c_prm,gradMu_prm,tenJ_prm,'o','MarkerSize',3);
xlabel('C'); ylabel('gradMu'); zlabel('J');
axis tight
subplot(1,2,2)
plot(c_prm,Mu_prm,'o','MarkerSize',3);
xlabel('C'); ylabel('Mu')
axis tight

% ASSIGNING THE DATA-SET TO THE GAUSS POINTS

% % Randomly assiging s\inS .
% gradMu_str = max(gradMu_prm)*rand(tne,ngp);
%   tenJ_str = max(tenJ_prm)*rand(tne,ngp);
%      c_str = max(c_prm)*rand(tne,ngp);
%     Mu_str = max(Mu_prm)*rand(tne,ngp);

% Constant s\inS .
gradMu_str = max(gradMu_prm)    *zeros(tne,ngp);
  tenJ_str = max(tenJ_prm)      *zeros(tne,ngp);
     c_str = min(c_prm)         *ones(tne,ngp);
    Mu_str = min(Mu_prm)        *ones(tne,ngp);
    
%% BOUNDARY CONDITION

% On chemical potential
pm = [1];                           % prescribed nodes
fm = setdiff(1:tnn,pm); fm = fm';   % free nodes

% On Lagrange multipliers
pl = pm;                            % prescribed nodes
fl = setdiff(1:tnn,pl); fl = fl';   % free nodes

% Initialize Mu, C and Lam
% Chemical potential field on the nodes
Mu  = min(Mu_prm)*ones(tnn,length(tItrt));
% Lagrange multiplier field on the nodes
Lam = zeros(tnn,length(tItrt));
% Concentration field on the Gauss integration points
C   = min(c_prm)*ones(tne,ngp,length(tItrt));

% Prescribing Mu in time.
% ----->
% % Sinusoidal Loading
% osc  = 1;
% a    = pi*osc*tItrt/T;
% Muosc = 1*sin(a);
% for n = 1 : length(tItrt)
%     Mu(1,n) = Muosc(n);
% end

% Constant Loading
Mu(1,:) = mu1;
% <-----

%% FECore
% In the following code the so-called stiffness and the mass matrices are assembled for once and for all. Even if the
% underlying data-set is made from a non-linear constitutive model. The data-driven solver does not require to reconstruct
% those matrices.
run('FECore_DataDriven.m');

Am =  barM+     barK;
Al =  -barM-dT^2*barK;

% P A T I T I O N I N G 
% ----->
Amfp = Am(fm,pm); Amff = (Am(fm,fm)); Fmf = barFm(fm,1); Fgf = barFg(fm,1);
Alfp = Al(fl,pl); Alff = (Al(fl,fl)); Fjf = barFj(fl,1); Fcf = barFc(fl,1); 

% Creating the list of Gauss points to plot concentration field
xgp=zeros(ngp,tne);
for kk = 1 : tne
    for ll = 1:ngp
        xgp(ll,kk) = N(ll,:)*x(egnn(kk,:));
    end
end

%% Data-Driven Solver

maxIter = 100;          % Maximum Iteration.
tol     = 1e-12;        % Tolerance and distance value. 
prvDis  = 1;            % Initializing previous distance value

for n = 2 : length(tItrt)   % Time stepping loop
    
    disp(['Time Step: ' num2str(n)])
    
    k = 1;                  % Iteration counter for data-driven solver 
    while k < maxIter       % Main data driven loop

        % Solve for Mu and Lam
        % ----->
        Mu (fm,n) = Amff \ ( Fmf + Fgf - Amfp*Mu (pm,n));
        Lam(fl,n) = Alff \ ( dT*Fjf - Fcf - Alfp*Lam (pl,n));
        % <-----
    
        % Initializing the force vector
        Dis = 0;
        barFg = zeros(tnn,1);
        barFj = zeros(tnn,1);
        barFc = zeros(tnn,1);
        barFm = zeros(tnn,1);
        
        % Element loop.
        for en = 1 : tne

            % Calling the global node numbering
            gnn = egnn(en,:);

            % Integration loop
            for gs = 1 : ngp
                
            % Jacobian Matrix
            Jcbn = B(gs,:)*x(egnn(en,:));
            % Calculating local states at (k) 
                % -------------------------------
                gradMu =                                 B(gs,:)/Jcbn * Mu (gnn,n) ;
                  tenJ = tenJ_str(en,gs) + dT*Dref*(B(gs,:)/Jcbn * Lam(gnn,n));
                  gpMu =                                 N(gs,:)      * Mu (gnn,n) ;
            C(en,gs,n) =     c_str(en,gs)     - 1/Lref * N(gs,:)      * Lam(gnn,n) ;
            
                %  DISTANCE
                Pi = dot((gpMu-Mu_prm)*1/Lref/nnb,(gpMu-Mu_prm),1)' + ...
                	 dot((C(en,gs,n)-c_prm)*Lref*nnb,(C(en,gs,n)-c_prm),1)' + ...
                     dot((gradMu-gradMu_prm) * Dref*nna,(gradMu-gradMu_prm),1)' + ...
                     dot((tenJ-tenJ_prm)*(1/Dref/nna),(tenJ-tenJ_prm),1)';

                 % FINDING THE POINT IN THE DATA-SET WHICH MINIMIZES THIS DISTANCE FUNCTION
                 [gpDis, indx] = min(Pi);
                 Dis = Dis + gpDis * glw(gs) * Jcbn;    % Integrating the distance globally
                 
                 % Assiging new internal-state 
                gradMu_str(en,gs) = gradMu_prm(indx);
                  tenJ_str(en,gs) =   tenJ_prm(indx);
                    Mu_str(en,gs) =     Mu_prm(indx);
                     c_str(en,gs) =      c_prm(indx);
                
                % Calculating the so-called force vectors
            	barFg(gnn,1) = barFg(gnn,1) + B(gs,:)'/Jcbn * Dref * gradMu_str(en,gs) * glw(gs) * Jcbn;
                barFm(gnn,1) = barFm(gnn,1) + N(gs,:)'    * 1/Lref *     Mu_str(en,gs) * glw(gs) * Jcbn;
                barFj(gnn,1) = barFj(gnn,1) + B(gs,:)'/Jcbn        *   tenJ_str(en,gs) * glw(gs) * Jcbn;
                barFc(gnn,1) = barFc(gnn,1) + N(gs,:)'    * (c_str(en,gs)-C(en,gs,n-1))  * glw(gs) * Jcbn;   
                 
            end
        end
        
        % Applying the boundary conditions on them
        Fmf = barFm(fm,1); Fgf = barFg(fm,1); Fjf = barFj(fl,1); Fcf = barFc(fl,1);
        
        % Post-processing
        disp(['    Iteration: ' num2str(k) '    Tolerance: ' num2str(abs(Dis-prvDis)) ]);
        
        figure(2)
        subplot(3,1,1);
        plot(x,Lam(:,n)); title('\lambda'); 
        ylim([-1 1]); xlim([min(x) max(x)])
        drawnow;
        figure(2)
        subplot(3,1,2);
        plot(x,Mu(:,n)); title('\mu');
        ylim([0 max(Mu_prm)]); xlim([min(x) max(x)])
        drawnow;
        figure(2)
        subplot(3,1,3);
        hold off
        for ll = 1 : ngp
            plot(xgp(ll,:),C(:,ll,n),'ko'); title('c at gauss point')
            hold on
        end
        ylim([min(c_prm) max(c_prm)]); xlim([min(x) max(x)])
        drawnow;
        
        % Check convergence
        % THERE IS SOME PROBLEM HERE BECAUSE THE INITIAL VALUE IS ZERO FURTHER IN THE SIMULATION.
        if abs(Dis-prvDis) < tol && k>5         
            break;
        else
            k = k + 1;
            prvDis = Dis;
        end
        
        
    end
    
    disp(' ');
    
end