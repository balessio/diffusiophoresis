function index_of_closest = ind_closest(val, array)
[~,index_of_closest] = min(abs(array-val));
end