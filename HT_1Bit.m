function [ zht ] = HT_1Bit( A,y,s )
% This function performs one iteration of hart thresholding 
% for sign measurements

%%% Backprojection
    zht = zeros(size(A,2),1);
      z = A' * y;
%%% Hard Thresholding
    [val,pos] = sort(abs(z),'descend');
    zht(pos(1:s),:) = z(pos(1:s),:);
end

