function [ dire_vec,dire_vec0, dire_vec1,G_rotated0 ] = edge_direct( G_rotated, G_ini,L )
% Calculate the direction for the L longest edges of graph "G_rotated"
% compute edge directions from acoustic ranging

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[n,m]=size(G_rotated);
dist_matrix=zeros(n);
dire_vec=zeros(L,m);
dire_vec0=zeros(L,m);
dire_vec1=zeros(L,m);
G_rotated0=G_rotated.*[-ones(n,1) ones(n,1)];%Inverted case

%compute the dist matrix of graph "G_rotated"
for i=1:n
    for j=1:i-1
        dist_matrix(i,j)=sqrt((G_rotated(i,:)-G_rotated(j,:))*(G_rotated(i,:)-G_rotated(j,:))');
    end
end

%calcute the direction vector
[~,Y]=sort(dist_matrix(:),'descend');

for i=1:L
    dire_vec(i,:)=G_rotated(mod(Y(i)-1,n)+1, :)-G_rotated(1+floor((Y(i)-1)/n), :);
    dire_vec0(i,:)=G_ini(mod(Y(i)-1,n)+1, :)-G_ini(1+floor((Y(i)-1)/n), :);
    dire_vec1(i,:)=G_rotated0(mod(Y(i)-1,n)+1, :)-G_rotated0(1+floor((Y(i)-1)/n), :);
end

end

