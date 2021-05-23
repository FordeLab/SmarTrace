function A = accumcells(cells, subfunc, field, valfun, fillval)
    subs = cell2mat(cellfun(subfunc, cells, 'UniformOutput', false));
    val = 1:length(cells);
    %fun = @(x) valfun(cells(x).{field});
    fun = @(x) valfun(cellfun(@(c) c.(field), cells(x)));
    A = accumarray(subs', val, [], fun, fillval);
