function remove_column = bt_remove_index(target_year_path,smos_father_path,sic_format)
    sic_files = dir([target_year_path,'*',sic_format]);
    target_freezep_range = length(sic_files);

    % 加2 是因为要将经纬度的那两列也计算进去
    remove_column = false(1,target_freezep_range+2);
    for j = 1:target_freezep_range
        sic_name = sic_files(j).name;
        
        pattern = '\d{8}'; % 匹配连续的8个数字
        sic_date = cell2mat(regexp(sic_name, pattern, 'match'));
        
        % 根据日期判断是否存在亮温文件
        smos_bt_path = [smos_father_path,'/daily_merge/'];
        smos_bt_name = [sic_date,'.txt'];
        if ~exist ([smos_bt_path,smos_bt_name],'file')
            remove_column(j+2) = true;
        end
    end 
end