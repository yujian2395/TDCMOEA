function DA = DiversityArchive(DA, Offspring, zmin, Ns, gen, maxGen)
    S  = [DA, Offspring];                  
    Obj = S.objs;                              
    Con = S.cons;                               
    CV  = sum(max(0, Con), 2);             
    NonDominated = DominationCal(S, 1,1e-6);
 
    S   = S(NonDominated);
    Obj = Obj(NonDominated, :);
    CV  = CV(NonDominated);

    [W, Ns] = UniformPoint(Ns, size(Obj, 2));   
    Angle_W_to_W = acos(1 - pdist2(W, W, 'cosine'));
    Angle_W_to_W(eye(Ns) == 1) = inf;
    h = mean(min(Angle_W_to_W));              
    if sum(CV <= 0) > 0
        zmax = max(Obj(CV <= 0, :), [], 1);
    else
        [~, index] = min(CV);
        zmax = Obj(index, :);
    end
    Obj_norm = (Obj - zmin) ./ (zmax - zmin + eps);
    cos_vals = 1 - pdist2(W, Obj_norm, 'cosine');
    cos_vals = max(cos_vals, -1);
    cos_vals = min(cos_vals, 1);
    Angle_S_to_W = sin(acos(cos_vals));       

    [w1, w2] = adaptive_weights(gen, maxGen, CV);

    DA = [];
 
    for i = 1 : Ns
        angles = Angle_S_to_W(i, :);
        neighbor = (angles <= h);
        
        if ~any(neighbor)
            [~, idx] = min(angles);
            selected = idx;
        else
            cand_indices = find(neighbor);
            cand_angles = angles(neighbor);
            cand_CV = CV(neighbor);
            
            [~, ~, angle_rank] = unique(cand_angles);
            [~, ~, cv_rank] = unique(cand_CV);
            
            combined_rank = w1 * cv_rank + w2 * angle_rank;
            [~, min_idx] = min(combined_rank);
            selected = cand_indices(min_idx);
        end
        
        DA = [DA, S(selected)];
    end
end

function [w1, w2] = adaptive_weights(gen, maxGen, CV)
    progress = gen / maxGen;  
    w1_base =0.1 + 0.2 * progress; 
    feasible_ratio = sum(CV <= 0) / length(CV); 
    adjust = 0.5 *  (1 - feasible_ratio);
    w1 = min(max(w1_base + adjust, 0.1), 1);
    w2 = 1 - w1;
end
