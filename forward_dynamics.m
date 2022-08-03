function [obe,groups,coordinates] = forward_dynamics(TI,ST,b,sigma,h,sigma_B,mean,cov)
%TI: length of a time interval
%ST: steps

temp = mvnrnd(mean,cov);
% Simulation starts
group = expmapping(temp(1:3));
coordinate = temp(4:6)';

obe_dim = 6;
obe = zeros(obe_dim,ST);
groups = zeros(3,3,ST);
coordinates = zeros(3,ST);
groups(:,:,1) = group;
coordinates(:,1) = coordinate;
obe(:,1) = 0;
for i=1:ST-1
    obe(:,i+1) = obe(:,i)+ h(i*TI,group,coordinate)*TI+sigma_B*sqrt(TI)*randn(obe_dim,1);
    
    group_next = group * expmapping(coordinate * TI);   
    coordinate_next = coordinate + b(i*TI,group,coordinate,zeros(3,1))*TI + sigma(i*TI,group,coordinate)*sqrt(TI)*randn(3,1);
    group = group_next;
    coordinate = coordinate_next;
    groups(:,:,i+1) = group;
    coordinates(:,i+1) = coordinate;
end

    









