function [obj_num] = pt_nms(ex_points, delta_dist)
    if isempty(ex_points)
        obj_num = 0;
    else
        x = ex_points(:, 1);
        y = ex_points(:, 2);
        v = ex_points(:, 3);
        [v, idxs] = sort(v);
        pick = [];
        while ~isempty(idxs)
            last = length(idxs);
            i_idxs = idxs(last);
            pick = [pick; i_idxs];
            suppress = [last];
            for pos = 1 : last - 1
                j_idxs = idxs(pos);
                dx = abs(x(i_idxs) - x(j_idxs));
                dy = abs(y(i_idxs) - y(j_idxs));
                if norm([dx, dy]) < delta_dist
                    suppress = [suppress; pos];
                end
            end
            idxs(suppress) = [];
        end
        obj_num = length(pick);
    end
end