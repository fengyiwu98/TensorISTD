function [neg_con_pix_num, pos_obj_num] = analyse_pts(gt_pts, detect_pts, area_list, radius)
    % gt_pts: GT�е�Ŀ���
    % detect_pts: ������Ŀ��λ��
    % area_list: ÿ��Ŀ�������б�
    % neg_con_pix_num: ��Ŀ���������ͨ��������
    % pos_obj_num: ��Ŀ������Χ�ڵ���ͨ�����
    % radius: �жϼ����GT��֮��ľ����Ƿ��㹻С
    if isscalar(area_list)
        area_list = [area_list];
    end
    dist = zeros(length(gt_pts) / 2, length(detect_pts) / 2);
    row = 0;
    for idx1 = 1 : 2 : length(gt_pts) % ��ʵ
        y_gt = gt_pts(idx1); % col
        x_gt = gt_pts(idx1 + 1); % col; % row
        row = row + 1;
        col = 0;
        for idx2 = 1 : 2 : length(detect_pts) % Ԥ��
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