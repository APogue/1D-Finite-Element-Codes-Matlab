% Returns Gauss Legendre integration points and weights for 2 and three point integration.
% 
% Created:       27 August, 2017
% Last Modified: 23 Nov, 2020
% Author: Abdullah Waseem

if ngp == 1
	
	glz(1,1) = 0;
	glw(1,1) = 2;
	
elseif ngp == 2
    
    glz(1,1) = -1/sqrt(3);
    glz(2,1) =  1/sqrt(3);
    
    glw(1,1) =  1;
    glw(2,1) =  1;

elseif ngp == 3
        
    glz(1,1) = -sqrt(3/5);
    glz(2,1) =  0;
    glz(3,1) =  sqrt(3/5);
    
    glw(1,1) =  5/9;
    glw(2,1) =  8/9;
    glw(3,1) =  5/9;
    
end

