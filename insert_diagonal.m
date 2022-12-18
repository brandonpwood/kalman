% be careful not to add anything to to offest (ie; inserting in the first
% element would be x = 0, y = 0)
function acc = insert_diagonal(target, x, y, value)
    s = size(value);
    acc = target;
    for i = [1:s(1)]
        for j = [1:s(2)]
            value(i, j);
            acc(x+i, y+j);
            acc(x+i, y+j) = value(i, j);
        end
    end
end