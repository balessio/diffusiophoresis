% factor tells you the fraction of the max to look for the x-loc of
function [width,xMax,xLeft,xRight,maxm] = width_nAvg(x,nAvg)
[maxm,indMax] = max(nAvg);
[~,indLeft] = min(abs(nAvg(1:indMax)-((maxm+nAvg(1))/2)));
[~,indRight] = min(abs(nAvg(indMax:end)-((maxm+nAvg(end))/2)));
indRight = indRight + indMax - 1;
xMax = x(indMax); xLeft = x(indLeft); xRight = x(indRight);
width = xRight - xLeft;
end