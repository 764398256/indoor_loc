function [ G_orient ] = orientation( dire_vec,dire_vec0,dire_vec1,G_rotated,G_rotated0 )
% Rotated G_rotated to find the optimal graph orientation that maximize the
% inner product summation between dire_vec and dire_vec0
% "G_orient" is the graph after orientation
% "G_rotated0": consider the invert case

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Construct inner product summation function
[m, ~]=size(G_rotated);
syms x x0;
g=[cos(x) -sin(x);sin(x) cos(x)];
g0=[cos(x0) -sin(x0);sin(x0) cos(x0)];
j=sum(sum(dire_vec*g.*dire_vec0));
j0=sum(sum(dire_vec1*g0.*dire_vec0));
h=inline(diff(sum(sum(dire_vec*g.*dire_vec0))));
h0=inline(diff(sum(sum(dire_vec1*g0.*dire_vec0))));
% k=inline(diff(sum(sum(dire_vec*g.*dire_vec0)),2));
% k0=inline(diff(sum(sum(dire_vec1*g0.*dire_vec0)),2));
%Find the optimal graph orientation
x=fzero(h,0);
x0=fzero(h0,0);
if abs(eval(j)) > abs(eval(j0))
    if eval(j)>0
        x=fzero(h,0);
    else
        x=fzero(h,0) + pi;
    end
    G_orient=(G_rotated - repmat(G_rotated(1,:), m, 1))*eval(g);%rotated around the 1st node
else
    if eval(j0)>0
        x0=fzero(h0,0);
    else
        x0=fzero(h0,0) + pi;
    end
    G_orient=(G_rotated0 - repmat(G_rotated(1,:), m, 1))*eval(g0);%rotated around the 1st node
    
end

