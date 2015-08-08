function [ A0 ] = find_components( dist_constr )
%%find the largest components of the graph
%return the ver number
dist_measured = dist_constr;
dist_measured(find(dist_constr == inf)) = 0;
[m,n]=size(dist_measured); % in this algorithm m=n

if m~=n
    disp('The distance matric is not squre, data error');
    return;
end

[A, B] = components(sparse(dist_measured));
A0 = find(A==1);
end

