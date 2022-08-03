function output = LieAlgebra(coordinates)
a = coordinates(1);
b = coordinates(2);
c = coordinates(3);
b1 = [0,0,0;0,0,-1;0,1,0];
b2 = [0,0,1;0,0,0;-1,0,0];
b3 = [0,-1,0;1,0,0;0,0,0];

output = a*b1+b*b2+c*b3;

