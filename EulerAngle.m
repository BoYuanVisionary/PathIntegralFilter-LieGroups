function angles = EulerAngle(SO3)
%This function returns the EulerAngle of a given ratation matrix.
%See the first proper Euler angles on the Wikipedia page 'Euler_Angles'
angles = zeros(3,1);
angles(1) = atan(SO3(3,1) / SO3(2,1));
angles(2) = atan(sqrt(1-SO3(1,1)^2) / SO3(1,1) );
angles(3) = atan(-SO3(1,3) / SO3(1,2));
end

