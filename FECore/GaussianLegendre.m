% Returns Gauss Legendre integration points and weights for 2 and three point integration.
% 
% Created:       27 August, 2017
% Last Modified: 11 March, 2018
% Author: Abdullah Waseem

if ngp == 1
	
	glz(1,1) = 0;
	glw(1,1) = 2;
	
elseif ngp == 2
    
    glz(1,1) = -0.5773502691896257;
    glz(2,1) =  0.5773502691896257;
    
    glw(1,1) =  1.0000000000000000;
    glw(2,1) =  1.0000000000000000;

elseif ngp == 3
        
    glz(1,1) = -0.7745966692414834;
    glz(2,1) =  0.0000000000000000;
    glz(3,1) =  0.7745966692414834;
    
    glw(1,1) =  0.5555555555555556;
    glw(2,1) =  0.8888888888888888;
    glw(3,1) =  0.5555555555555556;
    
end

