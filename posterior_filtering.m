function [final_weights, final_groups, final_coordinates,gamma] = posterior_filtering(TI,ST,SN,b,sigma,h,sigma_B,Y,mean,cov)
%TI: length of a time interval
%ST: steps
%SN: number of samples


u=@(t,g,xi) zeros(3,1);
ratio = 0.1; %H = ratio * ST
threshold = 0.5;
H = round(ratio * ST)+1; % rounding up, at least 2
final_groups = zeros(SN,ST,3,3);
final_coordinates = zeros(SN,ST,3);
final_weights = zeros(SN,ST);
%Initial values P(X_0 |Y_0) is just the prior distribution
for j=1:SN
    temp = mvnrnd(mean,cov);
    final_groups(j,1,:,:) = expmapping(temp(1:3));
    final_coordinates(j,1,:) = temp(4:6)';
end
final_weights(:,1) = 1/SN;

%Simulation starts The first stage is just to repeat smoothing 
for k=2:H-1
    groups = zeros(SN,k,3,3);
    coordinates = zeros(SN,k,3);
    for j=1:SN
        temp = mvnrnd(mean,cov);
        groups(j,1,:,:) = expmapping(temp(1:3));
        coordinates(j,1,:) = temp(4:6)';
        Initial_ratio = 1;
    end
    % Update weights and get new samples
    temp_weights = zeros(SN,1);
    for j=1:SN
        %Generate new samples
        for i=1:k-1
            %samples
            group = squeeze(groups(j,i,:,:));
            coordinate = squeeze(coordinates(j,i,:));
            groups(j,i+1,:,:) = group * expmapping(coordinate * TI);
            u_current = u(i*TI,group,coordinate);
            coordinates(j,i+1,:) = coordinate + b(i*TI,group,coordinate,u_current)*TI + sigma(i*TI,group,coordinate)*(u(i*TI,group,coordinate)*TI+sqrt(TI)*randn(3,1));
        end
        %weights
        for i=1:k-1
            group = squeeze(groups(j,i,:,:));
            coordinate = squeeze(coordinates(j,i,:));
            h_current = h(i*TI,group,coordinate);
            h_after = h(i*TI,squeeze(groups(j,i+1,:,:)),squeeze(coordinates(j,i+1,:)));
            GUdt1 = 1/(2*sigma_B^2)*norm(h_current,2)^2*TI + 1/sigma_B^2*dot(Y(:,i),h_after-h_current);
            GUdt2 = 1/2*norm(u_current,2)*TI;
            GUdt3 = dot(u_current,sqrt(TI)*rand(3,1));
            temp_weights(j)  = temp_weights(j)+ GUdt1+GUdt2+GUdt3;
        end
        temp_weights(j)  = temp_weights(j) - 1/sigma_B^2*dot(Y(:,H),h_after);
        temp_weights(j)  = exp(-temp_weights(j)) * Initial_ratio;
    end
    temp_weights = temp_weights / sum(temp_weights);
    final_groups(:,k,:,:) =  groups(:,k,:,:);
    final_coordinates(:,k,:,:) =  coordinates(:,k,:);
    final_weights(:,k) = temp_weights;
end

for k=H:H
    groups = zeros(SN,k,3,3);
    coordinates = zeros(SN,k,3);
    for j=1:SN
        temp = mvnrnd(mean,cov);
        groups(j,1,:,:) = expmapping(temp(1:3));
        coordinates(j,1,:) = temp(4:6)';
        Initial_ratio = 1;
    end
    % Update weights and get new samples
    temp_weights = zeros(SN,1);
    prior_weights = zeros(SN,1);
    for j=1:SN
        %Generate new samples
        for i=1:k-1
            %samples
            group = squeeze(groups(j,i,:,:));
            coordinate = squeeze(coordinates(j,i,:));
            groups(j,i+1,:,:) = group * expmapping(coordinate * TI);
            u_current = u(i*TI,group,coordinate);
            coordinates(j,i+1,:) = coordinate + b(i*TI,group,coordinate,u_current)*TI + sigma(i*TI,group,coordinate)*(u(i*TI,group,coordinate)*TI+sqrt(TI)*randn(3,1));
        end
        %weights
        for i=1:k-1
            group = squeeze(groups(j,i,:,:));
            coordinate = squeeze(coordinates(j,i,:));
            h_current = h(i*TI,group,coordinate);
            h_after = h(i*TI,squeeze(groups(j,i+1,:,:)),squeeze(coordinates(j,i+1,:)));
            GUdt1 = 1/(2*sigma_B^2)*norm(h_current,2)^2*TI + 1/sigma_B^2*dot(Y(:,i),h_after-h_current);
            GUdt2 = 1/2*norm(u_current,2)*TI;
            GUdt3 = dot(u_current,sqrt(TI)*rand(3,1));
            temp_weights(j)  = temp_weights(j)+ GUdt1+GUdt2+GUdt3;
            if  i==1
                prior_weights(j) = prior_weights(j) + GUdt1+GUdt2+GUdt3 - 1/sigma_B^2*dot(Y(:,2),h_after);
                prior_weights(j)  = exp(-prior_weights(j)) * Initial_ratio;
            end
        end
        temp_weights(j)  = temp_weights(j) - 1/sigma_B^2*dot(Y(:,k),h_after);
        temp_weights(j)  = exp(-temp_weights(j)) * Initial_ratio;
    end
    temp_weights = temp_weights / sum(temp_weights);
    final_groups(:,k,:,:) =  groups(:,k,:,:);
    final_coordinates(:,k,:,:) =  coordinates(:,k,:);
    final_weights(:,k) = temp_weights;
end

% The second stage is recursively performing filtering
% Initial prior infomation
prior_weights = prior_weights / sum(prior_weights);
prior_groups  = groups(:,2,:,:);
prior_coordinates = final_coordinates(:,2,:);
gamma = zeros(ST-H,1);
for k = (H+1):ST
    groups = zeros(SN,H,3,3);
    coordinates = zeros(SN,H,3);
    
    groups(:,1,:,:) = prior_groups;
    coordinates(:,1,:) = prior_coordinates;
    Initial_weights = prior_weights;

    % Update weights and get new samples
    temp_weights = zeros(SN,1);
    prior_weights = zeros(SN,1);
    for j=1:SN
        %Generate new samples
        for i=1:H-1
            %samples
            group = squeeze(groups(j,i,:,:));
            coordinate = squeeze(coordinates(j,i,:));
            groups(j,i+1,:,:) = group * expmapping(coordinate * TI);
            u_current = u(i*TI,group,coordinate);
            coordinates(j,i+1,:) = coordinate + b(i*TI,group,coordinate,u_current)*TI + sigma(i*TI,group,coordinate)*(u(i*TI,group,coordinate)*TI+sqrt(TI)*randn(3,1));
        end
        %weights
        for i=1:H-1
            group = squeeze(groups(j,i,:,:));
            coordinate = squeeze(coordinates(j,i,:));
            h_current = h(i*TI,group,coordinate);
            h_after = h(i*TI,squeeze(groups(j,i+1,:,:)),squeeze(coordinates(j,i+1,:)));
            GUdt1 = 1/(2*sigma_B^2)*norm(h_current,2)^2*TI + 1/sigma_B^2*dot(Y(:,k-H+i),h_after-h_current);
            GUdt2 = 1/2*norm(u_current,2)*TI;
            GUdt3 = dot(u_current,sqrt(TI)*rand(3,1));
            temp_weights(j)  = temp_weights(j)+ GUdt1+GUdt2+GUdt3;
            if  i==1
                prior_weights(j) = prior_weights(j) + GUdt1+GUdt2+GUdt3 - 1/sigma_B^2*dot(Y(:,2),h_after);
                prior_weights(j)  = exp(-prior_weights(j)) .* Initial_weights(j);
            end
        end
        temp_weights(j)  = temp_weights(j) - 1/sigma_B^2*dot(Y(:,k),h_after);
        temp_weights(j)  = exp(-temp_weights(j)) * Initial_weights(j);
    end
    temp_weights = temp_weights / sum(temp_weights);
    final_groups(:,k,:,:) =  groups(:,H,:,:);
    final_coordinates(:,k,:) =  coordinates(:,H,:);
    final_weights(:,k) = temp_weights;
    
    % Update prior information
    prior_groups  = groups(:,2,:,:);
    prior_coordinates = final_coordinates(:,2,:);
    prior_weights = prior_weights / sum(prior_weights);
    
    % Check resampling
    gamma(k-H) = 1/(SN * sum(final_weights(:,k).^2));
    if (gamma(k-H) < threshold)
        new_index = datasample(1:SN,SN,'Weights',final_weights(:,k));
        prior_groups = prior_groups(new_index,:,:);
        prior_groups = reshape(prior_groups,[SN,1,3,3]);
        prior_coordinates = prior_coordinates(new_index,:);
        prior_weights = prior_weights ./ final_weights(:,k);
        prior_weights = prior_weights / sum(prior_weights);
    end

    



    
end