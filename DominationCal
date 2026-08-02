function NonDominated = DominationCal(Population, Add, R)
    PopObj = Population.objs;
    PopObj = roundn(PopObj, -10);
    [N, M] = size(PopObj);

    PopCon = Population.cons;
    Cons   = sum(max(0, PopCon), 2);

   
    R_eff = R;
    if (Add == 0 || Add == 1) && R_eff > 0
        distMat = pdist2(PopObj, PopObj);
    end

    Dominated = false(1, N);
    neighborThresh = 10;   

    for i = 1 : N-1
        err = PopObj(i,:) - PopObj;          
        eq  = zeros(N, 1);
        max_err = max(err, [], 2);
        min_err = min(err, [], 2);

        for j = i+1 : N
            
            for k = 1 : M
                if err(j, k) ~= 0
                    break;
                end
                if k == M
                    eq(j) = 1;
                end
            end

            if Add == 0
            
                if eq(j) == 1
                    Dominated(j) = true;
                elseif min_err(j) >= 0
                    Dominated(i) = true;
                elseif max_err(j) <= 0
                    Dominated(j) = true;
                end
            else
         
                if eq(j) == 1
                    if Cons(i) < Cons(j)
                        Dominated(j) = true;
                    elseif Cons(i) > Cons(j)
                        Dominated(i) = true;
                    end
                elseif min_err(j) >= 0
                    if Cons(j) <= 0 || Cons(j) <= Cons(i)
                        Dominated(i) = true;
                    end
                elseif max_err(j) <= 0
                    if Cons(i) <= 0 || Cons(i) <= Cons(j)
                        Dominated(j) = true;
                    end
                end
            end
        end
    end

if (Add == 0 || Add == 1) && R_eff > 0

    overlapThresh = R_eff / 1000;   

    for i = 1 : N-1
        for j = i+1 : N

            if Dominated(i) || Dominated(j)
                continue;   
            end

            dij = distMat(i, j);

            if dij <= overlapThresh
                Dominated(j) = true;
                continue;
            end

            if dij <= R_eff

                localR = dij / 2;

                distThresh = dij / 5;

                neighborIdx_i = distMat(i,:) <= localR & distMat(i,:) > 0;
                neighborCount_i = sum(neighborIdx_i);

                if any(neighborIdx_i)
                    avgLocalDist_i = mean(distMat(i, neighborIdx_i));
                else
                    avgLocalDist_i = inf;
                end

                neighborIdx_j = distMat(j,:) <= localR & distMat(j,:) > 0;
                neighborCount_j = sum(neighborIdx_j);

                if any(neighborIdx_j)
                    avgLocalDist_j = mean(distMat(j, neighborIdx_j));
                else
                    avgLocalDist_j = inf;
                end

                if neighborCount_i + neighborThresh <= neighborCount_j
                    Dominated(j) = true;

                elseif neighborCount_i >= neighborCount_j + neighborThresh
                    Dominated(i) = true;

                elseif abs(neighborCount_i - neighborCount_j) < neighborThresh

                    if avgLocalDist_i >= avgLocalDist_j + distThresh
                        Dominated(j) = true;

                    elseif avgLocalDist_i + distThresh <= avgLocalDist_j
                        Dominated(i) = true;
                    end
                end
            end
        end
    end
end

    NonDominated = ~Dominated;
end
