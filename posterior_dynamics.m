function [weights, groups, coordinates] = posterior_dynamics(TI,ST,SN,b,sigma,h,sigma_B,Y,mean,cov)
%TI: length of a time interval
%ST: steps
%SN: number of samples
groups = zeros(SN,ST,3,3);
coordinates = zeros(SN,ST,3);
u=@(t,g,xi) zeros(3,1);
%Initial distribution

%Simulation starts
for j=1:SN   
    temp = mvnrnd(mean,cov);
    groups(j,1,:,:) = expmapping(temp(1:3));
    coordinates(j,1,:) = temp(4:6)';
    %Modify the initial_ratio if using another initial distribution
    Initial_ratio = 1;
end

% Update weights and get new samples
temp_weights = zeros(SN,1);
for j=1:SN
    %Genrate new samples
    for i=1:ST-1
        %samples
        group = squeeze(groups(j,i,:,:));
        coordinate = squeeze(coordinates(j,i,:));
        groups(j,i+1,:,:) = group * expmapping(coordinate * TI);
        u_current = u(i*TI,group,coordinate);
        coordinates(j,i+1,:) = coordinate + b(i*TI,group,coordinate,u_current)*TI + sigma(i*TI,group,coordinate)*sqrt(TI)*randn(3,1);
    end
    
    %Weights
    
    for i=1:ST-1
        group = squeeze(groups(j,i,:,:));
        coordinate = squeeze(coordinates(j,i,:));
        h_current = h(i*TI,group,coordinate);
        h_after = h(i*TI,squeeze(groups(j,i+1,:,:)),squeeze(coordinates(j,i+1,:)));
        GUdt1 = 1/(2*sigma_B^2)*norm(h_current,2)^2*TI + 1/sigma_B^2*dot(Y(:,i),h_after-h_current);
        GUdt2 = 1/2*norm(u_current,2)*TI;
        GUdt3 = dot(u_current,sqrt(TI)*rand(3,1));
        temp_weights(j)  = temp_weights(j)+ GUdt1+GUdt2+GUdt3;
    end
    temp_weights(j)  = temp_weights(j) - 1/sigma_B^2*dot(Y(:,ST),h_after);
    temp_weights(j)  = exp(-temp_weights(j)) * Initial_ratio;
end
temp_weights = temp_weights / sum(temp_weights);

weights = repmat(temp_weights,1,ST);


% function  [result] = iLQR(b,sigma,h,sigma_B,Y,g1,xi1,ST,TI)
%  % Initial guess
%  groups = zeros(ST,3,3);
%  algebras = zeros(ST,3);
%  groups(1,:,:) = g1;
%  algebras(1,:) = xi1;
%  for i=1:ST-1
%     group = squeeze(groups(i,:,:));
%     algebra = squeeze(algebras(i,:,:));
%     groups(i+1,:,:) =  squeeze(groups(i,:,:))  * expmapping(algebra);
%     algebras(i+1,:) = algebras(i,:) + b(i*TI,group,algebra)*TI + sigma(i*TI,group,algebra)*sqrt(TI)*randn(3,1);
%  end
%  
%  C_end = @(g,xi) -1/ sigma_B^2*dot(Y(:,ST),h(0,g,xi));
%  temp = squeeze(groups(ST,:,:));
%  x_T = [algebras(ST,:)';temp(:)];
%  [V_T,v_T] = finite_difference(C_end,x_T);
%  for i=ST-1:1
%      
%  end 
% end
% function [first_derivatives,second_derivatives]= finite_difference(target,x,flag==1)
% % Compute the first and second order derivatives
% % Assume x is a column vector
% % flag= 1 means  computing both derivatives
% step_size = 0.0001;
% dim = size(x);
% 
% % First-order
% first_derivatives = zeros(dim,1);
% for i = 1:dim
%     delta_x = zeros(dim,1);
%     delta_x = delta_x(i)+step_size;
%     first_derivatives(i) = (target(x+delta_x) - target(x-delta_x)) / (2*delta_x);
% end
% % Second_order
% second_derivatives = zeros(dim,dim);
% if isequal(flag,1)
%     for i = 1:dim
%         delta_x_i = zeros(dim,1);
%         delta_x_i = delta_x_i(i)+step_size;
%         for j=1:dim
%             delta_x_j = zeros(dim,1);
%             delta_x_j = delta_x_j(j)+step_size;
%             second_derivatives(i,j) = (target(x+delta_x_i+delta_x_j)+target(x)...
%                 -target(x+delta_x_j)-target(x+delta_x_i)) / (step_size^2);
%         end
%     end
% end
% end