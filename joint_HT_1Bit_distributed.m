function [ Zht ] = joint_HT_1Bit_distributed( A,y,n,m,L,s )
% This function jointly recovers L signals from separate linear sign measurements
% via one joint hard-thresholding step

%%% Backprojection
      Z = zeros(n,L);
    Zht = zeros(n,L);
    for k = 1:L
       Z(:,k) = A(m*(k-1)+1:k*m,:)' * y(:,k);
    end
%%% Hard thresholding
    row_norm = norms(Z,2,2);
    [val,pos] = sort(row_norm,'descend');
    Zht(pos(1:s),:) = Z(pos(1:s),:);
end

