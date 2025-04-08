function [neg_con_pix_num, pos_obj_num] = analyse_pts(gt_pts, detect_pts, area_list, radius)
    % gt_pts: GT中的目标点
    % detect_pts: 检测出的目标位置
    % area_list: 每个目标的面积列表
    % neg_con_pix_num: 非目标区域的连通域像素数
    % pos_obj_num: 在目标区域范围内的连通域个数
    % radius: 判断检测点和GT点之间的距离是否足够小
    if isscalar(area_list)
        area_list = [area_list];
    end
    dist = zeros(length(gt_pts) / 2, length(detect_pts) / 2);
    row = 0;
    for idx1 = 1 : 2 : length(gt_pts) % 真实
        y_gt = gt_pts(idx1); % col
        x_gt = gt_pts(idx1 + 1); % col; % row
        row = row + 1;
        col = 0;
        for idx2 = 1 : 2 : length(detect_pts) % 预测
            y_pr = detect_pts(idx2); % col
            x_pr = detect_pts(idx2 + 1); % col; % row
            col = col + 1;
            dist(row, col) = norm([y_gt - y_pr, x_gt - x_pr]);
        end
    end
    dist = dist < radius;
    dist2 = sum(dist, 1) > 0;
    pos_obj_num = sum(dist2);
    neg_con_pix_num = sum((1 - dist2) .* area_list);
end