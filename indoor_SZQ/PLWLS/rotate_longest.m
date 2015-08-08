function [ G_rotated, dist_constr] = rotate_longest( G_constr )
% Rotate "G_constr" such that the direction of the longest edge is parallel
% to the X axis

%  "G_totated" is the Graph after rotation

%*************************************************************************
[m, n]=size(G_constr);
dist_constr=zeros(m);

for i=1:m
    for j=i+1:m
        dist_constr(i,j)=sqrt((G_constr(i,:)-G_constr(j,:))*(G_constr(i,:)-G_constr(j,:))');
    end
end
dist_constr=dist_constr+dist_constr';

[rows,cols]=find(dist_constr==max(max(dist_constr)));
row=rows(1);
col=cols(1);
long_vec=G_constr(row,:)-G_constr(col,:);% the vector corresponding to the longest edge

rotate_M=[long_vec(1) -long_vec(2);long_vec(2) long_vec(1)]/sqrt(long_vec*long_vec');%Coordinate rotation

G_rotated=G_constr*rotate_M;

end

