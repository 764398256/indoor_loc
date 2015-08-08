%Construct distance matrix

m=length(Distance);

A=zeros(m-1,1);
B=zeros(m-1,1);
D=zeros(m-1,1);

for i=2:m
    temp=Devices{i,1};
    A(i-1)=temp(1)-48;
    B(i-1)=temp(3)-48;
    D(i-1)=Distance(i);
end
%Construct spare Matrix
S_Matrix=sparse(A,B,D,7,7);

D_Matrix=full(S_Matrix);
D_Matrix=D_Matrix+D_Matrix';
%clear m A B temp