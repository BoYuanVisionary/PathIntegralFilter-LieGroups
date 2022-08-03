function [relative_errors_algebras,relative_errors_groups_1,relative_errors_groups_2,relative_errors_groups_3,final_estimated_algrbras,final_estimated_groups] =evaluation(estimated_algebras,true_algebras,estimated_groups,true_groups,weights)

final_estimated_algrbras = squeeze(sum(estimated_algebras.*weights,1));
error_matrix = final_estimated_algrbras-true_algebras';
relative_errors_algebras = vecnorm(error_matrix,2,2) ./ vecnorm(true_algebras',2,2);

final_estimated_groups = squeeze(sum(estimated_groups.*weights,1));

relative_errors_groups_1 = zeros(size(final_estimated_groups,1),1);
relative_errors_groups_2 = zeros(size(final_estimated_groups,1),1);
relative_errors_groups_3 = zeros(size(final_estimated_groups,1),1);
% % Given by simply taking the average 
% for i = 1:size(final_estimated_groups,1)
%     %abs_error = norm(squeeze(pregroups(i,:,:))'* squeeze(true_groups(i,:,:)) - eye(3),1);
%     abs_error = trace(squeeze(final_estimated_groups(i,:,:))'* squeeze(true_groups(:,:,i)) - eye(3));
%     relative_errors_groups_1(i) = abs_error;
%     
% end
% % Given by log-exp mapping
% for j = 1:size(estimated_groups,2)
%     temp = zeros(3,3);
%     for i = 1:size(estimated_groups,1)
%         temp = temp + weights(i,j)*logm(squeeze(estimated_groups(i,j,:,:)));
%     end
%     temp = [temp(3,2),temp(1,3),temp(2,1)];
%     temp = expmapping(temp);
%     %abs_error = norm(temp'* squeeze(true_groups(j,:,:)) - eye(3),1);
%     abs_error = trace(temp'* squeeze(true_groups(:,:,j)) - eye(3));
%     relative_errors_groups_2(j) = abs_error;
% end
% 
% % Given by polar factors
% for i = 1:size(final_estimated_groups,1)
%     [U,~,V] = svd(squeeze(final_estimated_groups(i,:,:)));
% 
%     %abs_error = norm((U*V')'* squeeze(true_groups(i,:,:)) - eye(3),1);
%     abs_error = trace((U*V')'* squeeze(true_groups(:,:,i)) - eye(3));
%     relative_errors_groups_3(i) = abs_error;
% end




