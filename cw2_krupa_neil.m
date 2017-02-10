% ----- HYPER PARAMETERS ----- %
CrossoverProbability = 0.2;
MutationProbability = 1.0;
PopulationSize = 50;
Generations = 250;
% ----- HYPER PARAMETERS ----- %

% ----- GENERATE RANDOM FIRST POPULATION ----- %
Chromosome = [PopulationSize, 30];
for x = 1:PopulationSize
    for y = 1:3:30
        Chromosome(x,y) = randi([1, 4]);
        for z = y+1:y+3
            Chromosome(x,z) = randi([0, 9]);
        end
    end
    Chromosome(x,31) = 0;
end
% ----- GENERATE RANDOM FIRST POPULATION ----- %

Population = Chromosome;

%%%%% %%%%% %%%%% REPEAT BELOW FOR No.GENERATIONS %%%%% %%%%% %%%%%
for Generation=1:Generations
    % ----- MAKE EACH ANT RUN THE SIMULATION AND SAVE FITNESS LEVEL ----- %
    load muir_world.txt
    for x = 1:PopulationSize
        [Population(x,31), ~] = simulate_ant('muir_world.txt', Population(x,:));
    end
    % ----- MAKE EACH ANT RUN THE SIMULATION AND SAVE FITNESS LEVEL ----- %

    % ----- SORT THE ANTS FITNESSES ----- %
    Population = sortrows(Population, 31);
    % ----- SORT THE ANTS FITNESSES ----- %

    % ----- NORMALIZE THE FITNESSES ----- %
    NormalizedFitnesses = Population(:,31)/sum(Population(:,31));
    % ----- NORMALIZE THE FITNESSES ----- %

    % ----- CALCULATE THE CUMULATIVE SUM OF THE NORMALIZED FITNESSES ----- %
    SumNormalizedFitnesses = cumsum(NormalizedFitnesses);
    % ----- CALCULATE THE CUMULATIVE SUM OF THE NORMALIZED FITNESSES ----- %

    % ----- KEEP THE HIGHEST TWO, ELITE ANTS ----- %
    BestAnts = Population(PopulationSize-1:PopulationSize,:);
    % ----- KEEP THE HIGHEST TWO, ELITE ANTS ----- %

    % ----- CREATE NEW POPULATION AND ADD TWO HIGHEST ANTS ----- %
    for i=1:2
        NewPopulation(i, :) = BestAnts(i, :);
    end
    NewPopulationSize = 2;
    % ----- CREATE NEW POPULATION AND ADD TWO HIGHEST ANTS ----- %
    while NewPopulationSize < PopulationSize
            
    %{
    % ---------- USE RANK SELECTION ---------- %
    Population = sortrows(Population, 31);
    for i=1:PopulationSize-2
       Population(i,31) = i;
    end
    NormalizedFitnesses = Population(:,31)/sum(Population(:,31));
    SumNormalizedFitnesses = cumsum(NormalizedFitnesses);
    % ---------- USE RANK SELECTION ---------- %
        %}
    
       
    % ----- USE ROULETTE WHEEL TO SELECT TWO RANDOM CHROMOSOMES ----- %
        R1 = rand;
        R2 = rand;
        PopulationStop1 = find(SumNormalizedFitnesses>=R1, 1, 'first');
        TempChromosome_1 = Population(PopulationStop1, :);
        PopulationStop2 = find(SumNormalizedFitnesses>=R2, 1, 'first');
        TempChromosome_2 = Population(PopulationStop2, :);
    % ----- USE ROULETTE WHEEL TO SELECT TWO RANDOM CHROMOSOMES ----- %
    
    
   %{
    % ----- SINGLE POINT CROSSOVER WITH PROBABILITY CrossoverProbability ----- %
    if (rand < CrossoverProbability)
        RandomPoint = randi([1, 30]);
        TempChromosome_1_1(:, 1:RandomPoint) = TempChromosome_1(:, 1:RandomPoint);
        TempChromosome_2_2(:, 1:RandomPoint) = TempChromosome_2(:, 1:RandomPoint);
        TempChromosome_1_1(:, RandomPoint+1:31) = TempChromosome_2(:, RandomPoint+1:31);
        TempChromosome_2_2(:, RandomPoint+1:31) = TempChromosome_1(:, RandomPoint+1:31);
        TempChromosome_1 = TempChromosome_1_1;
        TempChromosome_2 = TempChromosome_2_2;
    end
    % ----- SINGLE POINT CROSSOVER WITH PROBABILITY CrossoverProbability ----- %
   %}
    
    
    % ----- TWO POINT CROSSOVER WITH PROBABILITY CrossoverProbability ----- %
    if (rand < CrossoverProbability)
        RandomPoint1 = randi([1, 15]);
        RandomPoint2 = randi([16, 30]);
        
        TempChromosome_1_1(:, 1:RandomPoint1) = TempChromosome_1(:, 1:RandomPoint1);
        TempChromosome_1_1(:, RandomPoint1+1:RandomPoint2) = TempChromosome_2(:, RandomPoint1+1:RandomPoint2);
        TempChromosome_1_1(:, RandomPoint2+1:31) = TempChromosome_1(:, RandomPoint2+1:31);
        
        TempChromosome_2_2(:, 1:RandomPoint1) = TempChromosome_2(:, 1:RandomPoint1);
        TempChromosome_2_2(:, RandomPoint1+1:RandomPoint2) = TempChromosome_1(:, RandomPoint1+1:RandomPoint2);
        TempChromosome_2_2(:, RandomPoint2+1:31) = TempChromosome_2(:, RandomPoint2+1:31);
        
        TempChromosome_1 = TempChromosome_1_1;
        TempChromosome_2 = TempChromosome_2_2;
    end
    % ----- TWO POINT CROSSOVER WITH PROBABILITY CrossoverProbability ----- %
    %}
    
    
    % ----- MUTATION WITH PROBABILITY MutationProbability ----- %
    if (rand < MutationProbability)
        MutationRand1 = randi([1, 30]);
        if (MutationRand1==1) || (mod((MutationRand1-1),3)==0)
            TempChromosome_1(:,MutationRand1) = randi([1,4]); 
        else
            TempChromosome_1(:,MutationRand1) = randi([0,9]);
        end
    end
    
    if (rand < MutationProbability)
        MutationRand2 = randi([1, 30]);
        if (MutationRand2==1) || (mod((MutationRand2-1),3)==0)
            TempChromosome_2(:,MutationRand2) = randi([1,4]); 
        else
            TempChromosome_2(:,MutationRand2) = randi([0,9]);
        end
        
        
    end
    % ----- MUTATION WITH PROBABILITY MutationProbability ----- %
    
    % ----- PUT INTO POPULATION IF ENOUGH SPACE ----- %
    if NewPopulationSize < PopulationSize
        NewPopulationSize = NewPopulationSize +1;
        NewPopulation(NewPopulationSize, :) = TempChromosome_1;
    end
    
    if NewPopulationSize < PopulationSize
        NewPopulationSize = NewPopulationSize +1;
        NewPopulation(NewPopulationSize, :) = TempChromosome_2;
    end
    % ----- PUT INTO POPULATION IF ENOUGH SPACE ----- %
    
    end
     % ----- DISPLAY THE HIGHEST FITNESSES AS GRAPH ----- %
    HighestFitness = Population(PopulationSize, :);
    FitnessPlot = HighestFitness(:,31);
    plot(Generation, FitnessPlot, '*'); hold on;
    axis([0 Generations 0 89]);
    % ----- DISPLAY THE HIGHEST FITNESSES AS GRAPH ----- %
    fprintf('Best Fitness Value in Generation %d is: %d \n', Generation, FitnessPlot);
    Population = NewPopulation;
end
% ----- SORT THE ANTS FITNESSES ----- %
Population = sortrows(Population, 31);
% ----- SORT THE ANTS FITNESSES ----- %
figure;
AntRoute = muir_world;
[fitness, trail] = simulate_ant('muir_world.txt', Population(PopulationSize,:));
for i = 1:size(trail)
    AntRoute(trail(i,1),trail(i,2)) = -1;
end
%%disp(fitness);
%%disp(trail);
%%disp(Population(PopulationSize,:));
imagesc(AntRoute);







