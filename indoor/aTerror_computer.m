function [amean,amedian,aninety] = aTerror_computer(err,tw)

err = reshape(err(:,tw:end),1,[]);
[h,stats] = cdfplot_ext(err);
close;
amean = stats.mean;
amedian = stats.median;
aninety = stats.ninety;

end