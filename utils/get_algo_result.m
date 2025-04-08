function get_algo_result(eval_algo_names, eval_data_names, img_types, ...
                         algo_base_path, data_base_path, res_base_path, time_path)

    for idx = 1 : length(eval_algo_names)
        algo_name = eval_algo_names{idx};
        addpath([algo_base_path, algo_name, '/']);
    end   
%     tab_idx = 1;
%     dataset = {};
%     method = {};
%     fps = {};

    for idx = 1 : length(eval_algo_names)

        algo_name = eval_algo_names{idx};
        algo = eval(algo_name);
        
        txt_path = [time_path, algo_name, '.txt'];
        fid = fopen(txt_path, 'a');
        
        for idx2 = 1 : length(eval_data_names)
 
            data_name = eval_data_names{idx2};
            
            disp(['#1 get_algo_result: ', algo_name, ' / ', data_name, '...']);
            
            seq_root_path = [data_base_path, data_name, '/'];
            img_dir = [];
%             res_cell = {};
            
            png_res_base_path = [res_base_path, data_name, '/', algo_name, '/'];
            
            if ~exist(png_res_base_path, 'dir')
                mkdir(png_res_base_path);
            end
            
            png_path = [res_base_path, data_name, '/', algo_name, '/*.png'];
            start_idx = length(dir(png_path)) + 1;

            for idx3 = 1 : length(img_types)
                type_cell = struct2cell(dir([seq_root_path, img_types{idx3}]));
                if isempty(type_cell)
                    continue
                end
                cell_name = type_cell(1, :);
                img_dir = [img_dir cell_name];
            end
            
            img_dir = sort(img_dir);

            if start_idx > length(img_dir)
                continue;
            end
            
            time = 0; 
            for idx3 = start_idx : length(img_dir)
                img_path = [seq_root_path, img_dir{idx3}];
                img = imread(img_path);
                if length(size(img)) > 2
                    img = rgb2gray(img);
                end
                
                tic %计算算法耗时 
                res = algo.process(im2double(img)).result; % input: 0~1
                time = time + toc;
                res_path = [png_res_base_path, num2str(idx3, '%04d'), '.png'];
                imwrite(res / (max(res(:)) + eps), res_path); % double型图像范围需要0~1才能顺利保存
            end

%             dataset{tab_idx} = data_name;
%             method{tab_idx} = algo_name;
           fps = 1.0 * (length(img_dir) - start_idx + 1) / time;
           times = 1.0 * time/(length(img_dir) - start_idx + 1) ;
%            fps = 1.0 * (length(img_dir)) / time;          

           time_name = datestr(now, 'yy-mm-dd_HH-MM-SS');
           fprintf(fid, '%s\n', [time_name, ':']);
           disp(['#1 ', algo_name, ' in ', data_name, ' run as ', num2str(times, '%.4f')]);
           fprintf(fid, '%s\n', [algo_name, ' in ', data_name, ' run as ', num2str(times, '%.4f'),'times.']);
%            tab_idx = tab_idx + 1;
        end
    end
    
    fclose(fid);
%     table_message = table(dataset, method, fps);
%     writetable(table_message, xlsx_path);
    
end

