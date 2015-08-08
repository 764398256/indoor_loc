function [ G_constr ] = MDS_MAP( dist_measured )
% Apply classical MDS to the distance matrix, retaining the firest 2 (or 3) 
% largest eigenvalues and eigenvectors to construct a 2-D (or 3-D) relative map.

% Where "G_constr" is the graph construct according to the pairwise ranging
% measurements; "dist_measured" repersent the pairwise acoustic ranging
% measurement

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[m,n]=size(dist_measured); % in this algorithm m=n

if m~=n
    disp('The distance matric is not squre, data error');
    return;
end
%execute floyd alg
for k = 1:n
    for i = 1:n
        for j = 1:n
            if dist_measured(i,k) + dist_measured(k,j) < dist_measured(i,j);   
                dist_measured(i,j) = dist_measured(i,k) + dist_measured(k,j);
            end
        end
    end
end
 
%Check connectivity of the graph
if ~isempty(find(dist_measured==Inf, 1))
    disp('The Graph is not connected, needs to be divided into connected subgraph');
    return;
end

%Apply clasical MDS to the distance matric(dist_measured)
H=eye(n)-ones(n)/n;
B=-0.5*H*dist_measured.^2*H;
[~, S, V]=svd(B);%Get B's eigenvector by using singular value decomposition
X=V*sqrt(S);
G_constr=X(:,1:2);% to construct a 2-Drelative map

end

