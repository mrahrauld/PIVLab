 n = 1;
 j = 11;
 A = reshape(result_conv_list{n}(:,:,j),[],1);
 [max1 ind1] = max(A)
 A(ind1) = -inf;
 max2 = max(A);
 max2/max1;
  A = reshape(result_conv_list2{n}(:,:,j),[],1);
 [max1 ind1] = max(A)
 A(ind1) = -inf;
 max2 = max(A);
 max2/max1;
%  
%  result_conv_sum = result_conv_list{1};
%  for i= 2:length(result_conv_list)
%      result_conv_sum =result_conv_sum + result_conv_list{i};
%  end
%  result_conv_sum = (result_conv_sum./max(max(result_conv_sum)))*255;
%  result_conv_sum = reshape(result_conv_sum(:,:,j),[],1);
%  [max1 ind1] = max(result_conv_sum)
%  