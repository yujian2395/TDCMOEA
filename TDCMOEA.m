classdef TDCMOEA < ALGORITHM
% <2026> <multi/many> <real/binary/permutation><constrained/none>
% Three-Stage Co-Evolutionary Framework Based on Neighborhood Density for Constrained Multi-objective Optimization Problems
% e1 --- 1 --- Type of environmental selection for forward exploration (1. SPEA2 2. NSGA-II 3. modified NSGA-III)
% e2 --- 1 --- Type of environmental selection for feasible exploitation (1. SPEA2 2. NSGA-II 3. modified NSGA-III)

    methods
        function main(Algorithm,Problem)
            [e1,e2] = Algorithm.ParameterSet(1,1);
            Ns = floor(Problem.N/3);
            Population = Problem.Initialization();
            maxGen = ceil(Problem.maxFE / Problem.N);
            gen = 0;
            FEA = Population;
            FA = [];
            DA = [];
            LPCA = [];
            Offspring = Population;
            zmin = min(Population.objs,[],1) - 1e-6;
            flag = 1;
            number_cons =  ConstraintGrouping(Population);
            number = number_cons{1};
            Mg = ceil(300/ length( number_cons));
            nextIdx = 2;
            while Algorithm.NotTerminated(FEA)
                gen = gen + 1;

                if gen/maxGen >= 0.35
                    flag = 2;
                    if mod(gen, Mg) == 0 && nextIdx <= length(number_cons)
                        number = number_cons{nextIdx};   
                        nextIdx = nextIdx + 1;
                    end
                end
                if gen/maxGen >= 0.7
                    flag = 3;
                end
                zmin = min(zmin, min(Offspring.objs, [], 1) - 1e-6);
                if flag == 1
                    FA = ForwardexplorationArchive(FA,Offspring, zmin, Ns, e1);
                    FEA = FeasibilityexploitationArchive(FEA, Offspring, Problem.N,e2); 
                    DA = DiversityArchive(DA, Offspring, zmin, Ns, gen, maxGen);
                    Pop3 = FEA(unidrnd(length(FEA), [1, floor(Problem.N/3)]));
                    Pop1 = FA;
                    Pop4 = DA;
                    MatingPool_Pop1_1 = randperm(length(Pop1));
                    MatingPool_Pop1_2 = randperm(length(Pop1));
                    MatingPool_Pop3_1 = randperm(length(Pop3));
                    MatingPool_Pop3_2 = randperm(length(Pop3));
                    MatingPool_Pop4_1 = randperm(length(Pop4));
                    MatingPool_Pop4_2 = randperm(length(Pop4));
                    if rand() < 0.5
                        Offspring1 = OperatorDE(Problem, Pop1, Pop1(MatingPool_Pop1_1), Pop1(MatingPool_Pop1_2),{1,0.5,1,1});
                        Offspring3 = OperatorDE(Problem, Pop3, Pop3(MatingPool_Pop3_1), Pop3(MatingPool_Pop3_2),{1,0.5,1,1});
                        Offspring4 = OperatorDE(Problem, Pop4, Pop4(MatingPool_Pop4_1), Pop4(MatingPool_Pop4_2),{1,0.5,1,1});
                    else
                        Offspring1 = OperatorGA(Problem, Pop1(MatingPool_Pop1_1), {1,20,1,1});
                        Offspring3 = OperatorGA(Problem, Pop3(MatingPool_Pop3_1), {1,20,1,1});
                        Offspring4 = OperatorDE(Problem, Pop4, Pop4(MatingPool_Pop4_1), Pop4(MatingPool_Pop4_2),{1,0.5,1,1});
                    end
                        Offspring = [ Offspring3,Offspring4,Offspring1];
                elseif  flag == 2
                    FEA = FeasibilityexploitationArchive(FEA, Offspring, Problem.N,e2); 
                    DA = DiversityArchive(DA, Offspring, zmin, Ns, gen, maxGen);
                    LPCA = LayeredProgressiveConstraintArchive(LPCA, Offspring, number , Ns );
                    Pop3 = FEA(unidrnd(length(FEA), [1, floor(Problem.N/3)]));
                    Pop2 = DA;
                    Pop4 = LPCA;
                    MatingPool_Pop2_1 = randperm(length(Pop2));
                    MatingPool_Pop2_2 = randperm(length(Pop2));
                    MatingPool_Pop3_1 = randperm(length(Pop3));
                    MatingPool_Pop3_2 = randperm(length(Pop3));
                    MatingPool_Pop4_1 = randperm(length(Pop4));
                    MatingPool_Pop4_2 = randperm(length(Pop4));
                    if rand() < 0.5
                        Offspring2 = OperatorDE(Problem, Pop2, Pop2(MatingPool_Pop2_1), Pop2(MatingPool_Pop2_2),{1,0.5,1,1});
                        Offspring3 = OperatorDE(Problem, Pop3, Pop3(MatingPool_Pop3_1), Pop3(MatingPool_Pop3_2),{1,0.5,1,1});
                        Offspring4 = OperatorDE(Problem, Pop4, Pop4(MatingPool_Pop4_1), Pop4(MatingPool_Pop4_2),{1,0.5,1,1});
                    else
                        Offspring2 = OperatorGA(Problem, Pop2(MatingPool_Pop2_1), {1,20,1,1});
                        Offspring3 = OperatorGA(Problem, Pop3(MatingPool_Pop3_1), {1,20,1,1});
                        Offspring4 = OperatorDE(Problem, Pop4, Pop4(MatingPool_Pop4_1), Pop4(MatingPool_Pop4_2),{1,0.5,1,1});
                    end

                    Offspring = [ Offspring2,Offspring3,Offspring4];
                else
                    DA = DiversityArchive(DA, Offspring, zmin, Problem.N/2, gen, maxGen);
                    FEA = FeasibilityexploitationArchive(FEA, Offspring, Problem.N,e2);
                    Pop1 = DA;
                    Pop2 = FEA(unidrnd(length(FEA), [1, floor(Problem.N/2)]));
                    MatingPool_Pop1_1 = randperm(length(Pop1));
                    MatingPool_Pop1_2 = randperm(length(Pop1));
                    MatingPool_Pop2_1 = randperm(length(Pop2));
                    MatingPool_Pop2_2 = randperm(length(Pop2));
                    if rand() < 0.5
                        Offspring1 = OperatorDE(Problem, Pop1, Pop1(MatingPool_Pop1_1), Pop1(MatingPool_Pop1_2),{1,0.5,1,1});
                        Offspring2 = OperatorDE(Problem, Pop2, Pop2(MatingPool_Pop2_1), Pop2(MatingPool_Pop2_2),{1,0.5,1,1});
                    else      
                        Offspring1 = OperatorGA(Problem, Pop1(MatingPool_Pop1_1), {1,20,1,1});
                        Offspring2 = OperatorGA(Problem, Pop2(MatingPool_Pop2_1), {1,20,1,1});
                    end

                    Offspring = [Offspring1,Offspring2];


                end
            end
        end
    end
end

