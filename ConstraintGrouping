function con_groups = ConstraintGrouping(Population)

    cons = Population.cons;
    [~, nCon] = size(cons);
    if nCon <= 1
        con_groups = {1:nCon};
        return;
    end
    V = cons > 0;
    J = zeros(nCon);
    for j = 1:nCon
        for k = j+1:nCon
            inter = sum(V(:,j) & V(:,k));
            union = sum(V(:,j) | V(:,k));
            if union == 0
                J(j,k) = 0;
            else
                J(j,k) = inter / union;
            end
            J(k,j) = J(j,k);
        end
    end
    J(eye(nCon)==1) = 1;
    distMat = 1 - J;
    distMat(eye(nCon)==1) = 0;
    Z = linkage(squareform(distMat), 'average');
    G = nCon;
    cluster_idx = cluster(Z, 'maxclust', G);
    con_groups = cell(1, G);
    for g = 1:G
        con_groups{g} = find(cluster_idx == g);
    end
end
