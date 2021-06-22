function [L,N,X,Y] = contour_exports(X,Y,Z)
M = contour(X,Y,Z);

levels = [];
nums_at_levels = [];
xData = {};
yData = {};

i=1;
size_M2 = size(M,2);
while i<size_M2
    levels = [levels;M(1,i)];
    nums_at_levels = [nums_at_levels;M(2,i)];
    xData{end+1} = M(1,i+1:i+M(2,i));
    yData{end+1} = M(2,i+1:i+M(2,i));
    i = i+M(2,i)+1;
end
L = levels; N = nums_at_levels; X = xData; Y = yData;
end
