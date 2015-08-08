function [mean, median, ninety] = Terror_computer(txcoor,tycoor,txalg,tyalg,numdev)

errorT = zeros(numdev,1);

for i=1:numdev
    errorT(i) =sqrt((txalg(i)-txcoor(i))*(txalg(i)-txcoor(i))+(tyalg(i)-tycoor(i))*(tyalg(i)-tycoor(i)));
end

[h,stats] = cdfplot_ext(errorT);
close;
mean = stats.mean;
median = stats.median;
ninety = stats.ninety;
end