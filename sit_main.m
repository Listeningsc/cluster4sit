function sit_main(Geo,sil,matrix,AP_Geo,divide)
    % 将脚本编写成函数进行统一
    
    %% 参数区 -----------------------------------------------------------------
    years  = [2010:2014,2016:2018];      
    method = 'AP4Kcenter';            % 当前方法

    multie_parent_path = '/scratch/lijiaxing/output/python_output/cluster/';
    bt_parent_path = '/project/lijiaxing/data_source/bt_data/L1C2epsg3413/';
    sic_father_path    = '/project/lijiaxing/data_source/SIC/';
    smos_father_path   = '/project/lijiaxing/data_source/bt_data/';
    
    if contains(divide,'3')
        divide_ratio = 0.3; 
    elseif contains(divide,'no')
        divide_ratio = 1;
    else
        num = regexp(divide, '\d+', 'match'); % 提取连续数字
        divide_ratio = 1/str2double(num); 
    end
    output_path = fullfile('/scratch/lijiaxing/output/matlab_output/cluster4sit',method,Geo,matrix,'Another_regular_pref',AP_Geo,sil,divide);
    if ~isfolder(output_path); mkdir(output_path); end

    %% 启动并行池 -------------------------------------------------------------
    if isempty(gcp('nocreate'))
        parpool(8);   % 全局池，只开一次
    end

    %% 主循环（按年份顺序执行） ---------------------------------------------
    for idx = 1:numel(years)
        target_year = years(idx);
        % multie_parent_path = '/project/lijiaxing/data_source/cluster/';
        if contains(method,'kmean')
            % 包含kmean的聚类方法
            cluster_parent_path = fullfile(multie_parent_path,method,num2str(target_year));
        elseif contains(method,'AP')
            % ap 聚类系点
            cluster_parent_path = fullfile(multie_parent_path,method,'grid_search',Geo,matrix,'Another_pref',AP_Geo,sil,num2str(target_year));
        else
            cluster_parent_path = fullfile(multie_parent_path,method,num2str(target_year));
        end

        % ---------- HDF/NC 格式判断 ----------
        if target_year == 2010
            sic_format = '.hdf';
        elseif target_year == 2011
            sic_format = '.nc';
        else
            sic_format = '.he5';
        end
        target_sic_path = fullfile(sic_father_path, sprintf('%d-%d',target_year,target_year+1),'/');
        remove_column   = bt_remove_index(target_sic_path, smos_father_path, sic_format);

        % ---------- 亮温序列 ----------
        bt_series_name = sprintf('%d_bt_series_reprocessed_interpolant2epsg3413_sep2April-mid.mat', target_year);
        bt_series_file = fullfile(bt_parent_path, bt_series_name);
        bt_series      = load(bt_series_file,'bt_series').bt_series;
        bt_series(:,remove_column) = [];
        bt_series_processed = sortrows(fillmissing(bt_series,'linear',2));

        % ---------- 系点 ----------
        if contains(method,'kmean')
            multie_name = sprintf('%d_kmedoids_best_centroid.csv', target_year);
            multie      = load(fullfile(cluster_parent_path,multie_name));
        elseif contains(method,'AP')
            multie_name = sprintf('%d_AP_centroids.csv', target_year);
            multie      = readmatrix(fullfile(cluster_parent_path,multie_name));
        else
            multie      = readCentroidCSVFiles(fullfile(cluster_parent_path,num2str(target_year)));
        end

        inx_mtie            = multie(:,8) > multie(:,9) - 0.02*multie(:,3);
        multie(inx_mtie,8)  = multie(inx_mtie,9) - 0.02*multie(inx_mtie,3);
        t_mtie         = multie;
        t_mtie(:,10)   = round((t_mtie(:,10)+t_mtie(:,11))/2,3);

        % ---------- 海冰厚度计算 ----------
        [~, bt_c] = size(bt_series_processed);
        sit_wm    = bt_series_processed;
        sit_wm_std    = bt_series_processed;
        sit_wm_3std    = bt_series_processed;
        for i = 3:bt_c          % 每日亮温列
            target_bt   = bt_series_processed(:,[1,2,i]);
            % —— 系点映射 ——
            inx_bt       = ismember(target_bt(:,1:2), t_mtie(:,1:2),'rows');
            mtie_bt   = target_bt(inx_bt,:);
            inx_tie   = ismember(t_mtie(:,1:2), mtie_bt(:,1:2),'rows');
            target_mtie = t_mtie(inx_tie,:);
            mtie_sit    = mtie_bt;

            inx_max  = mtie_bt(:,3) >= target_mtie(:,7) - target_mtie(:,8);
            inx_zero = mtie_bt(:,3) <= target_mtie(:,3) + target_mtie(:,4);
            inx_mid  = ~(inx_max | inx_zero);

            mtie_sit(inx_max,3)  = -1./target_mtie(inx_max,10) .* log(target_mtie(inx_max,8) ./ (target_mtie(inx_max,7) - target_mtie(inx_max,3)));
            mtie_sit(inx_zero,3) = 0;
            mtie_sit(inx_mid,3)  = -1./target_mtie(inx_mid,10) .* log((target_mtie(inx_mid,7) - mtie_bt(inx_mid,3)) ./ (target_mtie(inx_mid,7) - target_mtie(inx_mid,3)));

            % —— 非系点 ——
            target_bt(inx_bt,:) = [];
            
            sit_weightmean   = target_bt;
            sit_weightmean_std   = target_bt;
            sit_weightmean_3std   = target_bt;
            
            parfor j = 1:size(target_bt,1)   % ★ 内层并行
                if contains(method,'kmedoids')
                    [weight,mtie_used] = Kcenter3ties(target_bt(j,1:2), t_mtie);
                else
                    dist2bt = distance(target_bt(j,2),target_bt(j,1), t_mtie(:,2),t_mtie(:,1)) * 6371 * pi / 180;
                    dist2bt(dist2bt==0) = 1e-5;   % 避免 Inf
                    thres_dist      = quantile(dist2bt,divide_ratio);
                    inx_dist   = dist2bt < thres_dist;
                    weight  = 1./dist2bt(inx_dist);
                    mtie_used  = t_mtie(inx_dist,:);
                end

                st = zeros(numel(weight),1);
                for m = 1:numel(weight)
                    if target_bt(j,3) >= mtie_used(m,7) - mtie_used(m,8)
                        st(m) = -1/mtie_used(m,10)*log(mtie_used(m,8)/(mtie_used(m,7)-mtie_used(m,3)));
                    elseif target_bt(j,3) <= mtie_used(m,3) + mtie_used(m,4)
                        st(m) = 0;
                    else
                        st(m) = -1/mtie_used(m,10)*log((mtie_used(m,7)-target_bt(j,3))/(mtie_used(m,7)-mtie_used(m,3)));
                    end
                end     
                % 不剔除std
                sit_weightmean(j,3) = sum(st.*weight)/sum(weight);
                % 剔除单倍std
                thres_u_std = mean(st)+std(st);
                thres_d_std = mean(st)-std(st);
                inx_std = st>=thres_d_std & st<=thres_u_std;
                sit_weightmean_std(j,3) = sum(st(inx_std).*weight(inx_std))/sum(weight(inx_std));
                % 剔除3std
                thres_u_3std = mean(st)+3*std(st);
                thres_d_3std = mean(st)-3*std(st);
                inx_3std = st>=thres_d_3std & st<=thres_u_3std;
                sit_weightmean_3std(j,3) = sum(st(inx_3std).*weight(inx_3std))/sum(weight(inx_3std));   
            end
            target_sit     = sortrows([sit_weightmean; mtie_sit]);
            sit_wm(:,i)    = target_sit(:,3);
            
            target_sit_std     = sortrows([sit_weightmean_std; mtie_sit]);
            sit_wm_std(:,i)    = target_sit_std(:,3);
            
            target_sit_3std     = sortrows([sit_weightmean_3std; mtie_sit]);
            sit_wm_3std(:,i)    = target_sit_3std(:,3);
        end

        % ---------- 保存 ----------
        out_name = sprintf('%d_Fcycle_until_April2days_sit_based_%s.mat', target_year, method);
        save(fullfile(output_path,out_name),'sit_wm');
        
        out_name_std = sprintf('%d_Fcycle_until_April2days_sit_based_%s_std.mat', target_year, method);
        save(fullfile(output_path,out_name_std),'sit_wm_std');
        
        out_name_3std = sprintf('%d_Fcycle_until_April2days_sit_based_%s_3std.mat', target_year, method);
        save(fullfile(output_path,out_name_3std),'sit_wm_3std');
        
        fprintf('[Save] %d 年完成 → %s\n',target_year, out_name);
    end

    %% 关闭池（可选） ---------------------------------------------------------
    delete(gcp('nocreate'));
end

