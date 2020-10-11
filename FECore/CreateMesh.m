function [ L, lnn, nne, el, egnn, tnn, x ] = CreateMesh( elementtype, tne, xstart, xend  )
% Creates mesh of Line with 1D linear and quadratic elements. 
% Node numbering
%     Q1:         1*-------------*2
%     Q2:         1*------2------*3
% Important outputs:
%     egnn:       Element Global Node Numbering       (tne x lnn)  Connectivity Matrix
%     x:          Nodal Coordinates Vector
% 
% Created:       27 August, 2017
% Last Modified: 31 August, 2017
% Author: Abdullah Waseem

% Length of the domain: L
L = xend - xstart;

if elementtype == 'Q1'
    % Local node numbering
    lnn = 1 : 2;
    % Number of nodes in the element
    nne = size(lnn,2);
    % Length of the element.
    el = L / tne;
    % Total number of nodes in the mesh
    tnn = (tne*(size(lnn,2)-1))+1;
    egnn = zeros(tne, size(lnn,2));
    % Connectivity Matrix (egnn = Element Global Node Numbering)
    i = 1;
    for en = 1 : tne
        for I = lnn
            egnn(en,I) = i;
            i = i + 1;
        end
        i = i - 1;
    end
    % Nodal Coordinates
    x  = (xstart : el : xend)';
    
elseif elementtype == 'Q2'
    % Local node numbering
    lnn = 1 : 3;
    % Number of nodes in the element
    nne = size(lnn,2);
    % Length of the element.
    el = L / tne;
    % Total number of nodes in the mesh
    tnn = (tne*(size(lnn,2)-1))+1;
    egnn = zeros(tne, size(lnn,2));
    % Connectivity Matrix (egnn = Element Global Node Numbering)
    i = 1;
    for en = 1 : tne
        for I = lnn
            egnn(en, I) = i;
            i = i + 1;
        end
        i = i - 1;
    end
    % Nodal Coordinates
    x  = (xstart : el/2 : xend)';
    
end

