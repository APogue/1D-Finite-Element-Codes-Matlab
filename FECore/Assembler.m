function [ barK, barC, barF ] = Assembler( egnn, nne, tne, tnn, eK, eC, eF, astype )
% This function assembles the stiffness, mass matrices and force vector with either regular assembly
% process (prefered for small problems) or sparse assembly (prefered for large problems).
% 
% Author:           Abdullah Waseem       
% Created:          20-August-2017       
% Last Modified:    17-December-2017
switch astype
    % If assembler type is stndrd
    case 'standard'
    % Allocate space for K, M and F
    barK = zeros(tnn, tnn);        
    barC = zeros(tnn, tnn);        
    barF = zeros(tnn, 1);
    % Loop over all the elements
    for en = 1 : tne
        % Calling the global nodes in that element
        gn = egnn(en, :);
        % Loop over element nodes
        for I = 1 : nne
            for J = 1 : nne
                % Assembling K and M in a standard manner.      (gn = global node number)
                barK(gn(I), gn(J)) = barK(gn(I), gn(J)) + eK(I, J, en);
                barC(gn(I), gn(J)) = barC(gn(I), gn(J)) + eC(I, J, en);
            end
            % Assembling F in a standard manner.
            barF(gn(I), 1) = barF(gn(I), 1) + eF(I, 1, en);
        end
	end
	%
    % If the assembler type is sparse
    case 'sparse'
    % References: 
    % https://blogs.mathworks.com/loren/2007/03/01/creating-sparse-finite-element-matrices-in-matlab/#view_comments
    % RCV methos is used here.
    % Length of the C,R and V columns   len = (#nodes in X) x (#nodes in Y) x tne
    len = nne * nne * tne;
    % Initializing the location and storage vectors.
    Ridx = zeros(len, 1);       
    Cidx = zeros(len, 1);      
    VK = zeros(len, 1);         
    VC = zeros(len, 1);
    % Indexer
    idx = 1;
    % The force vector will be assembled in a standard manner -- Initializing barF
    barF = zeros(tnn, 1);
    % Loop over all the elements
    for en = 1 : tne
        % Calling the global nodes in that element
        gn = egnn(en, :);
        % Loop over the element nodes
        for I = 1 : nne
            for J = 1 : nne
                % Storing
                Ridx(idx, 1)   = gn(I);
                Cidx(idx, 1)   = gn(J);
                VK(idx, 1)  = eK(I, J, en);         % CRV adds up all the repeated indices
                VC(idx, 1)  = eC(I, J, en);         % CRV adds up all the repeated indices
                % Incrementing the indexer
                idx = idx + 1;
            end
            % Assembling F in a standard manner.
            barF(gn(I), 1) = barF(gn(I), 1) + eF(I, 1, en);
        end
    end
    % Sparse K and M matrices.
    barK = sparse(Ridx, Cidx, VK);        
    barC = sparse(Ridx, Cidx, VC);
end
% 
end