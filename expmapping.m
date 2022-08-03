function group = expmapping(coordinates)
%This function is to compute the exponetial mapping for SO(3). The method
%is given by the Rodrigues rotation formula.
%Input:a,b,c represent three coordinates. See details in https://www.cis.upenn.edu/~cis610/geombchap14.pdf Lemma 14.2.3 
a = coordinates(1);
b = coordinates(2);
c = coordinates(3);
theta = norm([a,b,c],2);
if isequal(theta,0)
    group = eye(3);
    return
end
A = LieAlgebra([a,b,c]);
B = [a;b;c].*[a,b,c;a,b,c;a,b,c];
group = cos(theta)*eye(3)+ sin(theta)/theta*A+(1-cos(theta))/theta^2*B;
%check if group is in SO(3)
error_threshold = 0.001;
if (abs(det(group)- 1)>error_threshold) || (norm((group*group')-eye(3))>error_threshold)
    message1 = 'The numerical error for exp mapping is too large for';
    message2 = [num2str(a),',',num2str(b),',',num2str(c)];
    warning([message1,message2]);
end
if isnan(group(1,1))
    warning('NaN warning');
end
    