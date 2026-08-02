function FA = ForwardexplorationArchive(FA, Offspring, zmin, Ns, Add)
    
    FA = [FA, Offspring];
    NonDominated = DominationCal(FA, 0, 1e-6);
    FA = FA(NonDominated);
    PopObj = FA.objs;
    [P, M] = size(PopObj);
    if P <= Ns
        FA = FA([1:P, unidrnd(P, [1, Ns-P])]);
        return;
    end
    
    [W, Ns] = UniformPoint(Ns, M);
    cosW = 1 - pdist2(W, W, 'cosine');
    cosW(eye(Ns) == 1) = 0;
    Angle_W = acos(cosW);
    Angle_W(eye(Ns) == 1) = inf;
    h = mean(min(Angle_W));
  
    zmax = max(PopObj, [], 1);
    PopObj_norm = (PopObj - zmin) ./ (zmax - zmin + 1e-10);
 
    Angle = acos(1 - pdist2(W, PopObj_norm, 'cosine')); 
    
    if Add == 2
        CrowdDis = CrowdingDistance(PopObj_norm);
    end

    if Add == 1 || Add == 3
        Distance = sqrt(sum(PopObj_norm.^2, 2)); 
    end
  
    selected = false(1, P);
    FA_new = [];
    for i = 1:Ns
       
        candidates = find(Angle(i, :) <= h & ~selected);
        if isempty(candidates)
        
            [~, idx] = min(Angle(i, ~selected));
            unselected_indices = find(~selected);
            best_idx = unselected_indices(idx);
        else
           
            switch Add
                case 1
                    weights = 1 ./ (Distance(candidates) + eps);
                case 2
                    weights = CrowdDis(candidates);
                case 3
                    fitness = Distance(candidates)' .* sin(Angle(i, candidates));
                    weights = 1 ./ (fitness + eps);
                otherwise
                    error('Unknown Add value');
            end
            probs = weights / sum(weights);
            idx = randsample(length(candidates), 1, true, probs);
            best_idx = candidates(idx);
        end
        FA_new = [FA_new, FA(best_idx)];
        selected(best_idx) = true;
    end
    FA = FA_new;
end
