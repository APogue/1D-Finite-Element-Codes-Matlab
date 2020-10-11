% Returns shape function and their derivatives for 1D FEM problem
% REFERENCE: http://www.me.metu.edu.tr/courses/me582/files/Handouts/Shape_Functions.pdf
% 
% Created:       27 August, 2017
% Last Modified: 11 March, 2018
% Author: Abdullah Waseem

if elementtype == 'Q1'
    % Initializing shape function and its gradient
    N  = zeros(ngp, 2);
    B = zeros(ngp, 2);
    % Gauss Legendre loop
    for gs = 1 : ngp
        % Evaluation of zeta = e
         e = glz(gs);
         % Shape functions:
         N(gs, 1) = 0.5 * (1 - e);
         N(gs, 2) = 0.5 * (1 + e);
         % Derivative of shape functions:
         B(gs, 1) = - 0.5;
         B(gs, 2) =   0.5;
    end
%     
elseif elementtype == 'Q2'
    % Initializing shape function and its gradient
    N = zeros(ngp, 3);
    B = zeros(ngp, 3);
    for gs = 1 : ngp
        % Evaluation of zeta = e
        e = glz(gs);
        % Shape functions:
        N(gs, 1) = 0.5 * e * (e-1);
        N(gs, 2) = -(e+1) * (e-1);
        N(gs, 3) = 0.5 * e * (e+1);
        % Derivative of shape functions:
        B(gs, 1) = 0.5 * (2 * e - 1);
        B(gs, 2) = - 2 * e;
        B(gs, 3) = 0.5 * (2 * e + 1);
    end
    
end

