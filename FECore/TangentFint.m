% INTERNAL FLUX VECTOR AND TANGENT STIFFNESS MATRIX. 
% In full NR scheme these two things can be calculated together for current iteration. 
% 
% Note: The internal flux vector and tangent matrix are very specific to PDE in hand and the constitutive model being used
% in this PDE. 
% 
% Author:           Abdullah Waseem       
% Created:          4-November-2018       
% Last Modified:    5-November-2018

% Initialize elemental Internal force vector and tangent matrix. 
Fint_e = zeros(nne,1,tne);
Te     = zeros(nne, nne, tne);

% Element loop
for en = 1 : tne
    % Gauss integration loop
	for gs = 1 : ngp
		
		% Jacobian Matrix
		Jcbn = B(gs,:)*x(egnn(en,:));
		
		% Flux Term: q = k0(1+Beta*U)*gradU
		q = k0(en)*( 1 + Beta(en)*N(gs,:)*un(egnn(en,:),1) ) * B(gs,:)/Jcbn*un(egnn(en,:),1);
		Fint_e(:,1,en) = Fint_e(:,1,en) + B(gs,:)'/Jcbn * q * glw(gs) * Jcbn;
		
		% Tangent Term: k0*dgradU + k0*Beta*U*dgradU + k0*Beta*gradU*dU
		K = k0(en)*B(gs,:)/Jcbn + k0(en)*Beta(en)*N(gs,:)*un(egnn(en,:),1)*B(gs,:)/Jcbn + k0(en)*Beta(en)*B(gs,:)/Jcbn*un(egnn(en,:),1)*N(gs,:);
		Te(:,:,en) = Te(:,:,en) + B(gs,:)'/Jcbn * K * glw(gs) * Jcbn;
		
	end
end

%Assembly
[ T, ~, FINT ] = Assembler( egnn, nne, tne, tnn, Te, Te, Fint_e, 'sparse' );
