function curves_drawer(fig_idx, eval_algo_names, eval_data_names, ...
    figure_base_path, mat_base_path, x_axis_ratio, FPR_thres)
% 12 color
    color_map = [ ...
%          55/255   126/255  184/255;  % 蓝
%          205/255  121/255  255/255;  % 紫
%          77/255   175/255  74/255;   % 绿
%          255/255  116/255  0/255;    % 橙
%          228/255  26/255   28/255;   % 红
%          255/255  217/255  47/255;   % 黄
%          140/255  86/255   75/255;   % 棕
%          255/255  108/255  145/255;  % 粉
%          153/255  153/255  153/255;  % 灰
%          102/255  194/255  165/255;  % 青绿
%          229/255  196/255  148/255;  % 卡其
%          100/255  190/255  240/255
         
         
       
         1 0 0;  % 红
         0 1 0;  % 绿
         0 0 1;  % 蓝
         205/255  121/255  255/255;  % 紫
         255/255  217/255  47/255;   % 黄
         0 1 1;  % 青
         255/255  116/255  0/255;    % 橙
         1 0 1;  % 品红
         255/255  217/255  47/255;   % 黄
         153/255  153/255  153/255;  % 灰
         229/255  196/255  148/255;  % 卡其
         255/255  108/255  145/255;  % 粉
         102/255  194/255  165/255;  % 青绿
         140/255  86/255   75/255;   % 棕
         

         
    ];
    LineType = {'-.',':','-.',':',':',':',':',':','-.',':','-.'};
    pot_map = {
%         '-h','','','','','','','','','',''
        '-*', '-d', '-o', '-x',  '-^','-v', '-^','->', '-<', '-p', '-s', '-+'
};
    
    for idx2 = 1 : length(eval_data_names) % different sequence
        data_name = eval_data_names{idx2};
        index_txt = [mat_base_path, '/index_results/', data_name, '.txt'];
        
        % roc curve
        curve_cell = {};
        leg_cell = {};
        AUC_list = [];
        AUC_list1 = [];
        AUC_list2 = [];
        % pr curve
        curve_cell2 = {};
        leg_cell2 = {};
        AP_list = [];
        F_list = [];
        WF_list = [];
        
        % roc2 curve
        curve_cell3 = {};
        leg_cell3 = {};
        AUC2_list = [];

        for idx = 1 : length(eval_algo_names) % different algorithm
            
            algo_name = eval_algo_names{idx};
            
            disp(['#3 roc_drawer: ', data_name, ' / ', algo_name, '...']);
            
            mat_path = [mat_base_path, '/curve_results/', data_name, '_', algo_name, '.mat'];
            load(mat_path);
            
            roc_relate = curve_collector_saver{1};
            pr_relate = curve_collector_saver{2};
            F_relate = curve_collector_saver{3};
            roc2_relate = curve_collector_saver{4};
            
            curve_tmp = roc_relate{1};
            % AUC_list = [AUC_list, -AUC]; % descend
            FPR_axis = curve_tmp(1, :);
            TPR_axis = curve_tmp(2, :);
            
            % FPR_axis filter
            table = [TPR_axis; FPR_axis];
            part = table(2, :) <= FPR_thres;
            part1 = table(:, part);
            pos1 = part1(:, 1);
            x1 = pos1(2);   
            y1 = pos1(1);
            mark = find(table(2, :) == x1 );
            if mark ~= 1
                pos2 = table(:, mark - 1);
                x2 = pos2(2);
                y2 = pos2(1);
                if ~(x1 <= FPR_thres && x2 >= FPR_thres)
                    error('postion of x1 x2 error!');
                end
                yMax = (y2 - y1) * (FPR_thres - x1) / (x2 - x1) + y1;
            else
                yMax = 1;
            end
            part2 =  [[yMax; FPR_thres], part1, [0;0]];
            
            for k=0:100
                t(k+1)=k/100;
            end
          
           % FPR_axis = part2(2, :); 
            %TPR_axis = part2(1, :);
            curve_cell{idx} = [FPR_axis; TPR_axis];
            AUC_list = [AUC_list, trapz(FPR_axis, TPR_axis)];
            AUC_list1 = [AUC_list1, trapz(t, TPR_axis)];
            AUC_list2 = [AUC_list2, trapz(t, FPR_axis)];
           
            AP = pr_relate{2};
            curve_cell2{idx} = pr_relate{1};
            AP_list = [AP_list, -AP]; % descend
            F_list = [F_list, F_relate{1}];
            WF_list = [WF_list, F_relate{2}];
            
            AUC2 = roc2_relate{2};
            curve_cell3{idx} = roc2_relate{1};
            AUC2_list = [AUC2_list, -AUC2]; % descend
        end
        
        OOP=-AUC_list+AUC_list1-AUC_list2;%OD
        TD=-AUC_list+AUC_list1;%TD
        BS=-AUC_list-AUC_list2;%BS
        TDBS=AUC_list1-AUC_list2;%BS
        
        SNPR=AUC_list1./AUC_list2;
        des_path=[figure_base_path ,data_name];
        if ~exist(des_path, 'dir')
            mkdir(des_path);
        end
        
        for idx = 1 : length(eval_algo_names) % different algorithm
            disp([data_name,'\',eval_algo_names{idx},':OOP=',num2str(OOP(idx)),',TD=',num2str(TD(idx)),...
                ',BS=',num2str(BS(idx)),',TDBS=',num2str(TDBS(idx)),',SNPR=',num2str(SNPR(idx)),',AUC_DF=',num2str(-AUC_list(idx)) ',AUC_Dt=',num2str(AUC_list1(idx)),...
                 ',AUC_Ft=',num2str(AUC_list2(idx))]);
        fid = fopen(index_txt, 'a');
        time_name = datestr(now, 'yy-mm-dd_HH-MM-SS');
        fprintf(fid, '%s\n', [time_name, ':']);
        fprintf(fid,'%s\n',[algo_name,':\',eval_algo_names{idx},' ',num2str(OOP(idx)),' ',num2str(TD(idx)),...
        ' ',num2str(BS(idx)),' ',num2str(TDBS(idx)),' ',num2str(SNPR(idx)),' ',num2str(-AUC_list(idx)) ,' ',num2str(AUC_list1(idx)),...
         ' ',num2str(AUC_list2(idx))]);
        fclose(fid);
        end  
        
        [mess, AUC_idx] = sort(AUC_list);
        [mess, AUC_idx1] = sort(AUC_list1,'descend');
        [mess, AUC_idx2] = sort(AUC_list2);
        [mess, AP_idx] = sort(AP_list);
        [mess, AUC2_idx] = sort(AUC2_list);

         fig=figure(1);     
         for idx7 = 1 : length(AUC_idx)
             %idx = idx7; %调参
             idx = AUC_idx(idx7); %算法比较
             AUC = -AUC_list(idx);
             FT = curve_cell{idx};
             FPR_axis = FT(1, :);
            TPR_axis = FT(2, :);
            plot3(FPR_axis,t,TPR_axis, 'Color', color_map(idx7,:), ...
               'LineStyle', LineType{idx7}, 'MarkerSize', 1, 'LineWidth', 3, 'MarkerSize', 3); 
           leg_cell{idx7} = [eval_algo_names{idx}];
           hold on;
         end
        grid on;
        view(-75,10);
        xlabel('FPR'); ylabel('\tau');zlabel('TPR');
        xlim([0 FPR_thres * x_axis_ratio]);
        ylim([0 1]);
        zlim([0 1]);
        %title('3D\_ROC');
        set(gca,'ZTick',(0:0.1:1));
        set(gca,'YTick',(0:0.2:1));
        set(gca,'XTick',(0:2*1e-5:1));
        legend( leg_cell, 'Location','east', 'Interpreter','latex');
        saveas(fig, [des_path,'/SOTA_1.png']);
        
        fig = figure(2);
        for idx7 = 1 : length(AUC_idx)
            %idx = idx7; %调参
            idx = AUC_idx(idx7); %算法比较
            AUC = AUC_list(idx);
            % AUC = AUC / FPR_thres; % AUC校正
            FT = curve_cell{idx};
            FPR_axis = FT(1, :);
            TPR_axis = FT(2, :);
            %pp = pot_map(idx7);
            plot(FPR_axis, TPR_axis, 'Color', color_map(idx7,:), ...
               'LineStyle', LineType{idx7}, 'MarkerSize', 1, 'LineWidth', 3, 'MarkerSize', 3);
           leg_cell{idx7} = ['[', num2str(-AUC, '%.5f'), ']', eval_algo_names{idx}];
            hold on;
        end
        grid on;
        xlabel('FPR');  
        ylabel('TPR');
        xlim([0 FPR_thres * x_axis_ratio]); % 将坐标范围限制在0~1
        ylim([0 1]);
        %title(['ROC of ', data_name], 'Interpreter','none');
        legend( leg_cell, 'Location','southeast', 'Interpreter','latex');
%         grid off;
        saveas(fig, [des_path, '/SOTA_2.png']);
        
        fig=figure(3);
        for idx7 = 1 : length(AUC_idx)
            %idx = idx7; %调参
            idx = AUC_idx2(idx7);
            AUC = AUC_list2(idx);
            % AUC = AUC / FPR_thres; % AUC校正
            FT = curve_cell{idx};
            FPR_axis = FT(1, :);
            TPR_axis = FT(2, :);
           plot(t, FPR_axis, 'Color', color_map(idx7,:), ...
               'LineStyle', LineType{idx7}, 'MarkerSize', 1, 'LineWidth', 3, 'MarkerSize', 3);
           leg_cell{idx7} = [ eval_algo_names{idx}];
           hold on;      
        end
        grid on;
        xlabel('\tau');  
        ylabel('FPR');
        xlim([0 1]);
        ylim([0 FPR_thres * x_axis_ratio]);
        %title('ROC of (\tau,FPR)');NorthEast
        legend(leg_cell, 'Location','southwest', 'Interpreter','latex');
        saveas(fig, [des_path, '/SOTA_3.png']);
           
        fig=figure(4);
        for idx7 = 1 : length(AUC_idx)
            %idx = idx7; %调参
            idx = AUC_idx1(idx7);
            AUC = AUC_list1(idx);
            % AUC = AUC / FPR_thres; % AUC校正
            FT = curve_cell{idx};
            FPR_axis = FT(1, :);
            TPR_axis = FT(2, :);
            plot(t, TPR_axis, 'Color', color_map(idx7,:), ...
               'LineStyle', LineType{idx7}, 'MarkerSize', 1, 'LineWidth', 3, 'MarkerSize', 3);
           leg_cell{idx7} = [ eval_algo_names{idx}];
           hold on;      
        end 
        grid on;
        xlabel('\tau');  
        ylabel('TPR');
        xlim([0 1]); % 将坐标范围限制在0~1
        ylim([0 1]);
        legend(leg_cell, 'Location','southwest', 'Interpreter','latex');    
        saveas(fig, [des_path, '/SOTA_4.png']);
        fig_idx = fig_idx + 1;
    end

end