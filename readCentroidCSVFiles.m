function data = readCentroidCSVFiles(folderPath)
    % 读取指定文件夹中包含 'centroid' 的 CSV 文件 (仅一行数据的时候)
    % 输入:
    %   folderPath - 包含CSV文件的文件夹路径
    % 输出:
    %   data - 包含读取数据的数组，每行为一个CSV文件的数据

    % 检查文件夹是否存在
    if ~isfolder(folderPath)
        error('指定的文件夹不存在: %s', folderPath);
    end

    % 获取文件夹中所有CSV文件
    files = dir(fullfile(folderPath, '*.csv'));

    % 过滤文件名中包含'centroid'的文件
    centroidFiles = files(arrayfun(@(x) contains(x.name, 'centroid'), files));

    % 初始化存储数据的数组
    data = [];

    % 遍历每个文件并读取
    for i = 1:length(centroidFiles)
        fileName = fullfile(folderPath, centroidFiles(i).name);
        fprintf('读取文件: %s\n', fileName);
        try
            % 读取CSV文件中的一行数据
            fileData = readmatrix(fileName);
            if size(fileData, 1) > 1
                warning('文件 %s 可能包含多行，仅保留第一行。', fileName);
                fileData = fileData(1, :);
            end
            % 合并数据（每个文件作为一行）
            data = [data; fileData];
        catch ME
            warning('文件 %s 读取失败: %s', fileName, ME.message);
        end
    end

    % 显示读取到的数据文件名
    if ~isempty(data)
        fprintf('成功读取 %d 个包含 "centroid" 的CSV文件。\n', size(data, 1));
    else
        warning('未找到包含 "centroid" 的CSV文件。');
    end
end
