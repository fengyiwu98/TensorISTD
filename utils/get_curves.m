function get_curves(eval_algo_names, eval_data_names, ...
                    thres_num, radius, res_base_path, ...
                    mat_base_path, txt_base_path, mask_base_path, preimg_type)
                
    % calculate
    obj_area = (radius * 2 + 1) * (radius * 2 + 1);
    for idx2 = 1 : length(eval_data_names) % different sequence
        
        data_name = eval_data_names{idx2};
        ann_path = [txt_base_path, data_name, '.txt'];
        %ann_path = [txt_base_path, '.txt'];%噪声序列测试
        
        % 构建mask
        data_mask_path = [mask_base_path, data_name, '/'];
        if ~exist(data_mask_path, 'dir')
            png_path = [res_base_path, data_name, '/', eval_algo_names{1}, '/',preimg_type];
            png_dir = struct2cell(dir(png_path));  
            name_list = {png_dir{1,:}};
            fig = fopen(ann_path);
            [mask_T_cell, obj_nums] = get_OBJ_mask(fig, length(name_list), ...
                res_base_path, data_name, eval_algo_names{1}, name_list, radius);
            fclose(fig);
        else
            [mask_T_cell, obj_nums] = get_OBJ_mask_pic(data_mask_path);
        end

        for idx = 1 : length(eval_algo_names) % different algorithm
 
            algo_name = eval_algo_names{idx};
            
            disp(['#2 get_roc_curves: ', data_name, ' / ', algo_name, '...']);

            data_algo_curve_mat = [mat_base_path, '/curve_results/'];
            
            if ~exist(data_algo_curve_mat, 'dir')
                mkdir(data_algo_curve_mat);
            end
            
            curve_mat_path = [data_algo_curve_mat, data_name, '_', algo_name, '.mat'];
            if exist(curve_mat_path, 'file')
                continue;
            end
            
            png_path = [res_base_path, data_name, '/', algo_name, '/', preimg_type];
            png_dir = struct2cell(dir(png_path));    
            name_list = {png_dir{1,:}};

            FPR_axis = zeros(1, thres_num + 1);
            TPR_axis = zeros(1, thres_num + 1);
            recall_axis = zeros(1, thres_num + 1);
            precision_axis = zeros(1, thres_num + 1);
            
            Pre_axis = zeros(1, thres_num + 1);
            Recall_axis = zeros(1, thres_num + 1);
            hitRate_axis = zeros(1, thres_num + 1);
            falseAlarm_axis = zeros(1, thres_num + 1);
            F_score = 0;
            WF_score = 0;
            tau=[0:0.01:1];
            % % Pixel values range from 0 to 1, so the threshold is
            for k = 0 : thres_num % - 1

                len =length(name_list);
                num_true = 0; num_false = 0; img_num = 0;
                FP = 0; TN = 0; TP = 0; FN = 0;
                Pre_S = 0; Recall_S = 0;
                hitRate_S = 0; falseAlarm_S = 0;

                for j = 1 : len

                    if obj_nums(j) == 0
                        continue;
                    end
                                               
                    % % Get an image, and obtain a logical image
                    preImg_path = [res_base_path, data_name, '/', algo_name, '/', name_list{j}];
                    preImg = im2double(imread(preImg_path));

                    if sum(preImg(:)) == 0
                        continue;
                    end
                    
                   % preImg_mask = (preImg >= tau(k+1));  % get threshold
                    preImg_mask = (preImg >= (k/(thres_num-1))); 
            
                    one_img_mask_Ts = mask_T_cell{j};
                    mask_T_all = zeros(size(one_img_mask_Ts{1}));
                    for idx9 = 1 : length(one_img_mask_Ts)
                        mask_T_all = logical(mask_T_all + one_img_mask_Ts{idx9});
                        tp1 = sum(sum(one_img_mask_Ts{idx9} .* preImg_mask)); % 真正
                        num_true = num_true + double(tp1 >= 1);
                    end
                    
                    if sum(mask_T_all(:)) == 0
                        continue;
                    end                    
                    
                    img_num = img_num + 1;
                    
                    % 目标区域内
                    tp = sum(sum(mask_T_all .* preImg_mask)); % 真正
                    fn = sum(sum(mask_T_all .* (~preImg_mask))); % 假负
                    TP = TP + tp;
                    FN = FN + fn;

                    % 背景区域内
                    mask_F = ~mask_T_all;
                    tn = sum(sum(mask_F .* (~preImg_mask))); % 真负
                    fp = sum(sum(mask_F .* preImg_mask)); % 假正
                    TN = TN + tn;
                    FP = FP + fp;  
                    num_false = num_false + fp;
                    
                    % pr-plus
                    if k == 10 % F-measure的计算和阈值无关所以计算一次即可
                        [PreFtem, RecallFtem, FmeasureF] = calculate_Fmeasure(preImg, mask_T_all);
                        [Weighted_FmeasureF] = calculate_Weighted_Fmeasure(preImg, mask_T_all);
                        F_score = F_score + FmeasureF;
                        WF_score = WF_score + Weighted_FmeasureF;
                    end
                    
                    [Pre, Recall, hitRate, falseAlarm] = Pre_Recall_hitRate(preImg_mask, mask_T_all);
                    Pre_S = Pre_S + Pre; 
                    Recall_S = Recall_S + Recall;
                    hitRate_S = hitRate_S + hitRate; 
                    falseAlarm_S = falseAlarm_S + falseAlarm;

                end
                
                [m1,m2] = size(preImg);
                
                FPR_axis(k + 1) = num_false/(img_num*(m1*m2-obj_area)); % 方法1: 1.0 * FP / (FP + TN + eps); 方法2 /num_false/(img_num*(m1*m2-obj_area));
                TPR_axis(k + 1) = num_true / sum(obj_nums(1:len));% 1.0 * TP / (TP + FN + eps);  方法二num_true / sum(obj_nums(1:len));
                recall_axis(k + 1) = 1.0 * TP / (TP + FN + eps);
                precision_axis(k + 1) = 1.0 * (TP + eps) / (TP + FP + eps);
                
                Pre_axis(k + 1) = Pre_S / img_num;
                Recall_axis(k + 1) = Recall_S / img_num;
                hitRate_axis(k + 1) = hitRate_S / img_num;
                falseAlarm_axis(k + 1) = falseAlarm_S / img_num;
                
                if k == 10 % F-measure的计算和阈值无关所以计算一次即可
                    F_score = F_score / img_num;
                    WF_score = WF_score / img_num;
                end
                
            end

            auc = -trapz(FPR_axis, TPR_axis);
            ap = -trapz(Recall_axis(1:end-1), Pre_axis(1:end-1));
            auc2 = -trapz(falseAlarm_axis, hitRate_axis);
            curve_collector_saver = {{[FPR_axis; TPR_axis], auc}, ...
                {[Recall_axis(1:end-1); Pre_axis(1:end-1)], ap}, ...
                {F_score, WF_score} ...
                {[falseAlarm_axis; hitRate_axis], auc2}}; % (data_name, algo_name)
            save(curve_mat_path, 'curve_collector_saver');
        end
    end  
end

function [mask_T_cell, obj_nums] = get_OBJ_mask(fig, len, res_base_path, data_name, ...
    algo_name, name_list, radius)

    mask_T_cell = {};
    obj_nums = [];
    
    for idx = 1 : len
        preImg_path = [res_base_path, data_name, '/', algo_name, '/', name_list{idx}];
        preImg = im2double(imread(preImg_path));
        img_sz = size(preImg);
        gt_line = fgetl(fig);
        if gt_line == -1
            continue;
        end
        if strcmp(gt_line, 'None')
            obj_nums = [obj_nums, 0];
            mask_T_cell{idx} = {-1};
            continue;
        end
        if strcmp(gt_line, '0 		0')
            obj_nums = [obj_nums, 0];
            mask_T_cell{idx} = {-1};
            continue;
        end
        gt_bboxes = strsplit(gt_line);
        gt_bboxes(cellfun(@isempty, gt_bboxes)) = []; % clear empty cell
        % 获取目标区域（正样本）掩模图
%         mask_T = zeros(img_sz);
        obj_num = 0;
        mask_Ts = {};
        idxx = 1;
        for idx3 = 1 : 2 : length(gt_bboxes) % all boxes
            y = round(str2double(gt_bboxes{idx3})); % col
            x = round(str2double(gt_bboxes{idx3 + 1})); % row
            
            if x <= radius || y <= radius || x >= img_sz(1) - radius || y >= img_sz(2) - radius
                continue;
            end
            obj_num = obj_num + 1;
            mask = zeros(img_sz);
            mask(x - radius : x + radius, y - radius : y + radius) = 1; % 如果在边界会改变原图尺寸
%             mask_T = logical(mask_T + logical(mask));
            mask_Ts{idxx} = logical(mask);
            idxx = idxx + 1;
        end
%         mask_T_cell{idx} = mask_T;
        mask_T_cell{idx} = mask_Ts;
        obj_nums = [obj_nums, obj_num];
    end
end

function [mask_T_cell, obj_nums] = get_OBJ_mask_pic(data_mask_path)

    mask_path = [data_mask_path, '/*.png'];
    mask_dir = struct2cell(dir(mask_path));  
    name_list = {mask_dir{1,:}};
    
    mask_T_cell = {};
    obj_nums = [];
    
    for idx = 1 : length(name_list)
        mask_T = {};
        mask_path2 = [data_mask_path, '/', name_list{idx}];
        mask = imread(mask_path2); % source image
        if length(size(mask)) > 2
            mask = rgb2gray(mask);
        end
        labeled_mask = bwlabel(logical(mask));
        obj_num = max(labeled_mask(:));
        obj_nums = [obj_nums, obj_num];
        for idx2 = 1 : obj_num
            e_mask = zeros(size(mask));
            e_mask(labeled_mask == idx2) = 1;
            mask_T{idx2} = logical(e_mask);
        end
        mask_T_cell{idx} = mask_T;
    end
end

function [Pre, Recall, hitRate, falseAlarm] = Pre_Recall_hitRate(testMap, gtMap)
    %testMap：输入图像与阈值比较后的逻辑矩阵
    %gtMap：真值逻辑矩阵
    %输出：
    neg_gtMap = ~gtMap; %取相反数
    neg_testMap = ~testMap;

    hitCount = sum(sum(testMap.*gtMap));%二值分割后的图像中占真值中的元素（1）的个数（交）
    trueAvoidCount = sum(sum(neg_testMap.*neg_gtMap));%既不属于真值也不属于二值后的图像的元素的个数（1）
    missCount = sum(sum(testMap.*neg_gtMap));%二值分割后图像中错分的1的个数
    falseAvoidCount = sum(sum(neg_testMap.*gtMap));%二值后，真值中没有没检测到的个数

    if hitCount==0
        Pre = 0;
        Recall = 0;
    else
        Pre = hitCount/(hitCount + missCount );
        Recall = hitCount/(hitCount + falseAvoidCount);
    end

    falseAlarm = 1 - trueAvoidCount / (eps+trueAvoidCount + missCount);
    hitRate = hitCount / (eps+ hitCount + falseAvoidCount);
end

function [PreFtem, RecallFtem, FmeasureF] = calculate_Fmeasure(sMap, gtMap)

    sumLabel = 2 * mean(sMap(:));
    
    if sumLabel > 1
        sumLabel = 1;
    end

    Label3 = zeros(size(gtMap));
    Label3(sMap >= sumLabel) = 1;

    NumRec = length(find(Label3==1));
    LabelAnd = Label3 & gtMap;
    NumAnd = length(find(LabelAnd==1));
    num_obj = sum(sum(gtMap));

    if NumAnd == 0
        PreFtem = 0;
        RecallFtem = 0;
        FmeasureF = 0;
    else
        PreFtem = NumAnd/NumRec;
        RecallFtem = NumAnd/num_obj;
        FmeasureF = ((1.3* PreFtem * RecallFtem ) / (0.3 * PreFtem + RecallFtem)); % beta^2=0.3
    end
end

function [Q] = calculate_Weighted_Fmeasure(FG, GT)
% WFb Compute the Weighted F-beta measure (as proposed in "How to Evaluate
% Foreground Maps?" [Margolin et. al - CVPR'14])
% Usage:
% Q = FbW(FG,GT)
% Input:
%   FG - Binary/Non binary foreground map with values in the range [0 1]. Type: double.
%   GT - Binary ground truth. Type: logical.
% Output:
%   Q - The Weighted F-beta score

    %Check input
    if (~isa( FG, 'double' ))
        error('FG should be of type: double');
    end
    if ((max(FG(:))>1) || min(FG(:))<0)
        error('FG should be in the range of [0 1]');
    end
    if (~islogical(GT))
        error('GT should be of type: logical');
    end

    dGT = double(GT); %Use double for computations.


    E = abs(FG-dGT);
    % [Ef, Et, Er] = deal(abs(FG-GT));

    [Dst,IDXT] = bwdist(dGT);
    %Pixel dependency
    K = fspecial('gaussian',7,5);
    Et = E;
    Et(~GT)=Et(IDXT(~GT)); %To deal correctly with the edges of the foreground region
    EA = imfilter(Et,K);
    MIN_E_EA = E;
    MIN_E_EA(GT & EA<E) = EA(GT & EA<E);
    %Pixel importance
    B = ones(size(GT));
    B(~GT) = 2-1*exp(log(1-0.5)/5.*Dst(~GT));
    Ew = MIN_E_EA.*B;

    TPw = sum(dGT(:)) - sum(sum(Ew(GT))); 
    FPw = sum(sum(Ew(~GT)));

    R = 1- mean2(Ew(GT)); %Weighed Recall
    P = TPw./(eps+TPw+FPw); %Weighted Precision

    Q = (2)*(R*P)./(eps+R+P); %Beta=1;
    % Q = (1+Beta^2)*(R*P)./(eps+R+(Beta.*P));
end
