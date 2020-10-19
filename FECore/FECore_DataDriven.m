%% 1D FEM CORE FOR TRANSIENT DATA-DRIVEN DIFFUSION MECHANICS

% Initializing
barFg = zeros(tnn,1);
barFj = zeros(tnn,1);
barFc = zeros(tnn,1);
barFm = zeros(tnn,1);
Ke = zeros(nne,nne,tne);
Me = zeros(nne,nne,tne);
Fe = zeros(nne,1,tne);

% Element loop
for en = 1 : tne
    
    gnn = egnn(en,:);
    
    % Gauss integration loop
	for gs = 1 : ngp
		
		% Jacobian Matrix
		Jcbn = B(gs,:)*x(egnn(en,:));
        
		% Element Stiffness Matrix
		Ke(:,:,en) = Ke(:,:,en) + B(gs,:)'/Jcbn * Dref * B(gs,:)/Jcbn * glw(gs) * Jcbn;
		
		% Element Mass Matrix
		Me(:,:,en) = Me(:,:,en) + N(gs,:)' * 1/Lref * N(gs,:) * glw(gs) * Jcbn;

	end
end
% Assemble barK, barC and barF
[ barK, barM, ~ ] = Assembler( egnn, nne, tne, tnn, Ke, Me, Fe, 'sparse' );