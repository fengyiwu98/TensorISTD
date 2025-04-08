function measure_calculator(eval_algo_names, eval_data_names, ...
                                           data_base_path, res_base_path, mat_base_path,...
                                           txt_base_path, img_types, preimg_type)
    % calculate
    rt = 2;     % radius of target
    rb = 10;     % radius of background
    
    
    for idx2 = 1 : length(eval_data_names) % different sequence
        data_name = eval_data_names{idx2};
        
        index_txt = [mat_base_path, 'index_results/', data_name, '.txt'];
        fid = fopen(index_txt, 'a');
        
        for idx = 1 : length(eval_algo_names) % different algorithm
            algo_name = eval_algo_names{idx};
            
            png_path = [res_base_path, data_name, '/', algo_name, '/', preimg_type];
            png_dir = struct2cell(dir(png_path));    
            name_list = {png_dir{1,:}};
            len = 90;%length(name_list);
            len = length(name_list);

            ann_path = [txt_base_path, data_name, '.txt'];
            %ann_path = [txt_base_path, '.txt'];%噪声序列测试
            fig = fopen(ann_path);
            
            seq_root_path = [data_base_path, data_name, '/'];

            % Calculate the sum of S and C
            sin = [];    sout = [];
            cin = [];    cout = [];
            BSR_array = [];

            % define mask to extract background
            mask_back = ones(2 * rb + 1, 2 * rb + 1);
            mask_back(rb - rt + 1 : rb + rt + 1, rb - rt + 1 : rb + rt + 1) = 0;
            mask_back = logical(mask_back);
            
            for idx3 = 1 : len
                
                % image name
                number = num2str(idx3, '%04d');
%                 number = num2str(idx3);
                
                for idx31 = 1 : length(img_types)
                    img_type = img_types(idx31);
                    img_type = img_type{1};
                    src_path = strcat(seq_root_path, number, img_type(2:end));
                    if ~exist(src_path, 'file')
                        continue;
                    end
                    img = imread(src_path); % source image
                    if length(size(img)) > 2
                        img = rgb2gray(img);
                    end
                    img = im2double(img);
                    break;
                end

              %  assert(exist('img', 'var') == true); % img得正确读取

                preImg_path = [res_base_path, data_name, '/', algo_name, '/', name_list{idx3}];
                preImg = im2double(imread(preImg_path));
                
                if sum(preImg(:)) == 0
                    continue;
                end
                
                % pad to avoid target in edge, now pass this step
                pad_radius = rb + 1;
                inImg = padarray(img, [pad_radius, pad_radius], 'replicate');
                outImg = padarray(preImg, [pad_radius, pad_radius], 'replicate');
                
                gt_line = fgetl(fig);
                
                if gt_line == -1
                    continue;
                end

                if strcmp(gt_line, 'None') % 图中无目标
                    continue;
                end
                if strcmp(gt_line, '0 		0') % 图中无目标
                    continue;
                end
                
                gt_bboxes = strsplit(gt_line);
                gt_bboxes(cellfun(@isempty, gt_bboxes)) = [];
                
                for idx31 = 1 : 2 : length(gt_bboxes) % all boxes
                    
                    try
                        y = round(str2double(gt_bboxes{idx31})) + pad_radius; % col
                        x = round(str2double(gt_bboxes{idx31 + 1})) + pad_radius; % row
                        % We define a matrix with a radius to be target and others are background
                        % get the background patch    
                        patch_in = inImg(x - rb : x + rb, y - rb : y + rb);
                        patch_out = outImg(x - rb : x + rb, y - rb : y + rb);
                    catch err
                        disp(['The ', num2str(idx3), 'th frame of ', data_name, ' has A edge target !']);
                        continue;
                    end

%                     max_inT = inImg(x, y);
%                     max_outT = outImg(x, y);
                    max_inT = max(max(patch_in(~mask_back)));
                    max_outT = max(max(patch_out(~mask_back)));

                    mean_inB = mean(patch_in(mask_back));
                    mean_outB = mean(patch_out(mask_back));

                    sin = [sin, abs(max_inT - mean_inB)];
                    sout = [sout, abs(max_outT - mean_outB)];

                    % calculate the standard deviation and assign to cin/cout
                    cin = [cin, std(patch_in(mask_back))];
                    cout = [cout, std(patch_out(mask_back))];
                    
                    BSR_array = [BSR_array, log(mean_inB / mean_outB)];
                end
            end
            
            fclose(fig);

%             SCR_Gain_array = (sout ./ cout) ./ (sin ./ cin);
            SCR_Gain_array = (sout ./ cout) ./ (sin ./ cin);
            BSF_array = cin ./ cout;
            CG_array = sout ./ sin;
            
            SCR_Gain_array(abs(SCR_Gain_array) == inf) = [];
            BSF_array(abs(BSF_array) == inf) = [];
            CG_array(abs(CG_array) == inf) = [];
            BSR_array(abs(BSR_array) == inf) = [];
            
            % here should be a check to NaN
            SCR_Gain = nanmean(SCR_Gain_array);
            BSF = nanmean(BSF_array);
            CG = nanmean(CG_array);
            BSR = nanmean(BSR_array);
            
            if isnan(SCR_Gain)
                SCR_Gain = inf;
            end
            if isnan(BSF)
                BSF = inf;
            end
            if isnan(CG)
                CG = inf;
            end
            if isnan(BSR)
                BSR = inf;
            end
            
            disp_str = ['#4 measure_calculator ', data_name, ' \ ', algo_name, ': SCR Gain = ', ...
                num2str(SCR_Gain,'%.4f'), ', BSF = ', num2str(BSF, '%.4f'),...
                ', CG = ', num2str(CG, '%.4f'), ', BSR = ', num2str(BSR, '%.4f')];
            
            disp(disp_str);
            fprintf(fid, '%s\n', disp_str(22:end));
            
%             dataset_name{tab_idx} = data_name;
%             method_name{tab_idx} = algo_name;
%             SCRG_val{tab_idx} = SCR_Gain;
%             BSF_val{tab_idx} = BSF;
%             tab_idx = tab_idx + 1;

        end
    end
    
    fclose(fid);
    
%     table_message = table(dataset_name, method_name, SCRG_val, BSF_val);
%     writetable(table_message, xlsx_path);
    
end