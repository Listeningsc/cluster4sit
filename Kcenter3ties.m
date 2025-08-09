function [weights,sel_3ties] = Kcenter3ties(bt_point, t_mtie)
    % weightedInterpolationKmeanKmedoids  对单个亮温点用 kmean_kmedoids 聚类结果插值
    %
    % 输入：
    %   bt_point : [lon, lat, bt_val]
    %   t_mtie   : N×M 聚类中心参数矩阵，第 1 列 lon，第 2 列 lat，
    %              第 3 列 h0，第 4 列 …，第 7 列 h1，第 8 列 h2，第 10 列 t0
    %
    % 输出：
    %   sit_val  : 加权插值后的海冰厚度

    lon_bt = bt_point(1);
    lat_bt = bt_point(2);

    % 1) 计算到所有系点的球面距离（km）
    dist2bt = distance(lat_bt, lon_bt, t_mtie(:,2), t_mtie(:,1))* 6371 * pi/180;
    dist2bt(dist2bt==0) = 1e-5;

    % 2) 找最近的点 A
    [~, sortIdx] = sort(dist2bt);
    idxA = sortIdx(1);

    % 3) 处理经度差，取 −180..+180
    lonA = t_mtie(idxA,1);
    all_lons = t_mtie(:,1);
    delta_lon = mod(all_lons - lonA + 180, 360) - 180;

    % 4) 东侧最近
    east_idxs = find(delta_lon>0);
    if isempty(east_idxs)
        [~, idxE] = min(delta_lon);
    else
        [~, rel] = min(delta_lon(east_idxs));
        idxE = east_idxs(rel);
    end

    % 5) 西侧最近
    west_idxs = find(delta_lon<0);
    if isempty(west_idxs)
        [~, idxW] = max(delta_lon);
    else
        [~, rel] = min(abs(delta_lon(west_idxs)));
        idxW = west_idxs(rel);
    end

    % 取三点
    sel_idx   = [idxA; idxE; idxW];
    sel_3ties  = t_mtie(sel_idx,:);

    d_sel = distance(lat_bt, lon_bt, sel_3ties(:,2), sel_3ties(:,1)) * 6371 * pi/180;
    d_sel(d_sel==0) = 1e-5;
    w        = 1 ./ d_sel;
    weights  = w / sum(w);
end
