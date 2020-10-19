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
		
        barFg(gnn,1) = barFg(gnn,1) + B(gs,:)'/Jcbn * Dref * gradMu_str(en,gs) * glw(gs) * Jcbn;
        barFm(gnn,1) = barFm(gnn,1) + N(gs,:)'    * 1/Lref *     Mu_str(en,gs) * glw(gs) * Jcbn;
        barFj(gnn,1) = barFj(gnn,1) + B(gs,:)'/Jcbn        *   tenJ_str(en,gs) * glw(gs) * Jcbn;
		barFc(gnn,1) = barFc(gnn,1) + N(gs,:)'    * (c_str(en,gs)-C(en,gs,1))  * glw(gs) * Jcbn;
        
		% Element Stiffness Matrix
		Ke(:,:,en) = Ke(:,:,en) + B(gs,:)'/Jcbn * Dref * B(gs,:)/Jcbn * glw(gs) * Jcbn;
		
		% Element Mass Matrix
		Me(:,:,en) = Me(:,:,en) + N(gs,:)' * 1/Lref * N(gs,:) * glw(gs) * Jcbn;

	end
end
% Assemble barK, barC and barF
[ barK, barM, ~ ] = Assembler( egnn, nne, tne, tnn, Ke, Me, Fe, 'sparse' );