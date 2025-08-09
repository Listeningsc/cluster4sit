function M = reconstructMatrix(N, indx)
%RECONSTRUCTMATRIX  根据已删列矩阵 N 与逻辑索引 indx 复原原矩阵 M
%
%   M = reconstructMatrix(N, indx)
%
%   输入
%   ----- 
%   N    : m×k   矩阵，已删除 indx==true 的列
%   indx : 1×n   或 n×1 逻辑向量，长度 n = 原始列数，
%                 true 表示该列在原矩阵中被删除 / 置为 NaN
%
%   输出
%   -----
%   M    : m×n   复原后的矩阵，indx==true 的列填充 NaN，其余列等于 N
%
%   例子
%   ----
%   M = reconstructMatrix(N, indx);

    % ---------- 输入检查 ----------
    if ~islogical(indx)
        error('indx 必须是逻辑向量 (logical array)。');
    end
    
    indx = indx(:).';                 % 转为行向量
    [m, k] = size(N);                 % N 的行列
    n      = numel(indx);             % 原始列数
    
    if sum(~indx) ~= k
        error(['indx 中 false 的数量 (%d) 与 N 的列数 (%d) 不一致。'], ...
               sum(~indx), k);
    end
    
    % ---------- 构造并填充 ----------
    M          = NaN(m, n, 'like', N);  % 与 N 同类 (double / single…) 的 NaN
    M(:,~indx) = N;                     % 仅在 false 位填充数据
end
