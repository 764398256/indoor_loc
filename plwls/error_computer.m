function [errDT,errT,error] =error_computer(xcoor,ycoor,xalg,yalg,numdev,numT)

errDT = zeros(numdev,numT);
errT = zeros(1,numT);
for i=1:numT
    for j=1:numdev
        dis =sqrt((xalg(j,i)-xcoor(j,i))*(xalg(j,i)-xcoor(j,i))+(yalg(j,i)-ycoor(j,i))*(yalg(j,i)-ycoor(j,i)));
        errDT(j,i) = dis;
    end
    errT(i) = mean(errDT(:,i));
end
error = mean(errT(:));
end