function Fpatch = bridge_rings_idx(inner_idx, outer_idx, V, orientation)
%BRIDGE_RINGS_IDX Zipper triangulation between two index loops.
%   Fpatch = BRIDGE_RINGS_IDX(inner_idx, outer_idx, V, orientation)
%   triangulates between two vertex index loops with robust handling of
%   winding order, degeneracies, and greedy diagonal selection.

    inner_idx = inner_idx(:)';
    outer_idx = outer_idx(:)';
    Ni = numel(inner_idx);
    No = numel(outer_idx);
    if Ni < 2 || No < 2
        Fpatch = zeros(0, 3);
        return;
    end

    center_point = mean(V(inner_idx, :), 1);
    w = orientation(:)' / norm(orientation);
    vecs_from_center = bsxfun(@minus, V(inner_idx, :), center_point);
    [~, max_dist_idx] = max(sum(vecs_from_center.^2, 2));
    ref_vec = V(inner_idx(max_dist_idx), :) - center_point;
    u = ref_vec - dot(ref_vec, w) * w;
    u = u / norm(u);
    v = cross(w, u);

    inner_vecs = bsxfun(@minus, V(inner_idx, :), center_point);
    angles_inner = atan2(inner_vecs * v', inner_vecs * u');
    [~, sort_order_inner] = sort(angles_inner);
    inner_idx = inner_idx(sort_order_inner);

    if No > 2
        sorted_outer_idx = zeros(1, No);
        dists_from_center_sq = sum(bsxfun(@minus, V(outer_idx, :), center_point).^2, 2);
        [~, start_pos] = max(dists_from_center_sq);
        sorted_outer_idx(1) = outer_idx(start_pos);
        remaining_indices = 1:No;
        remaining_indices(start_pos) = [];
        for k = 2:No
            last_added_v_idx = sorted_outer_idx(k - 1);
            remaining_v_indices = outer_idx(remaining_indices);
            coords_last_added = V(last_added_v_idx, :);
            coords_remaining = V(remaining_v_indices, :);
            distances_sq = sum(bsxfun(@minus, coords_remaining, coords_last_added).^2, 2);
            [~, closest_rem_pos] = min(distances_sq);
            sorted_outer_idx(k) = remaining_v_indices(closest_rem_pos);
            remaining_indices(closest_rem_pos) = [];
        end
        outer_idx = sorted_outer_idx;
    end

    v_inner_start = V(inner_idx(1), :);
    distances_sq = sum(bsxfun(@minus, V(outer_idx, :), v_inner_start).^2, 2);
    [~, best_start_j_pos] = min(distances_sq);
    outer_idx = circshift(outer_idx, [0, 1 - best_start_j_pos]);

    V_outer_on_plane = V(outer_idx, :) * [u', v'];
    signed_area = sum(V_outer_on_plane(1:end-1, 1) .* V_outer_on_plane(2:end, 2)) ...
                - sum(V_outer_on_plane(2:end, 1) .* V_outer_on_plane(1:end-1, 2));
    signed_area = signed_area + (V_outer_on_plane(end, 1) * V_outer_on_plane(1, 2) ...
                               - V_outer_on_plane(1, 1) * V_outer_on_plane(end, 2));

    if signed_area < 0
        outer_idx = [outer_idx(1), fliplr(outer_idx(2:end))];
    end

    Ftmp = zeros(Ni + No, 3);
    c = 0;
    i = 1; j = 1; ai = 0; aj = 0;

    while ai < Ni || aj < No
        curr_i_idx = inner_idx(i); curr_j_idx = outer_idx(j);
        inext = mod(i, Ni) + 1; jnext = mod(j, No) + 1;
        next_i_idx = inner_idx(inext); next_j_idx = outer_idx(jnext);
        if ai == Ni
            new_face = [curr_i_idx, next_j_idx, curr_j_idx];
            j = jnext; aj = aj + 1;
        elseif aj == No
            new_face = [curr_i_idx, next_i_idx, curr_j_idx];
            i = inext; ai = ai + 1;
        else
            d1 = sum((V(curr_j_idx, :) - V(next_i_idx, :)).^2);
            d2 = sum((V(curr_i_idx, :) - V(next_j_idx, :)).^2);

            if d1 < d2
                new_face = [curr_i_idx, next_i_idx, curr_j_idx];
                i = inext; ai = ai + 1;
            else
                new_face = [curr_i_idx, next_j_idx, curr_j_idx];
                j = jnext; aj = aj + 1;
            end
        end

        c = c + 1;

        p1 = V(new_face(1), :);
        p2 = V(new_face(2), :);
        p3 = V(new_face(3), :);
        face_normal = cross(p2 - p1, p3 - p1);

        if dot(face_normal, orientation) < 0
            Ftmp(c, :) = [new_face(1), new_face(3), new_face(2)];
        else
            Ftmp(c, :) = new_face;
        end
    end
    Fpatch = Ftmp(1:c, :);
end
