% Konstantinos Letros 8851
% Optimization Techniques
% The Project - Parameters Estimation
% Parameters Estimation using Genetic Algorithm

%% Clean the screen

clc
clear
close all;
format long;

tic
%% Function to be estimated

% Intput Vector: u
f = @(u) sin(u(1));
U=-2:0.05:2;
plot(U,arrayfun(f,U))
title("Original Function")
pause(0.01);

%% Initialization

% Number of Chromosomes in every Population
chromeNum = 50;

% Number of Gaussians in every Chromosome
gaussiansNum = 15;

% Number of generations until termination
generationsNum = 500;

% Initialize Population of Chromosomes
% population(gaussian_i,parameter_j,chromosome_k)
population = initPopulation(gaussiansNum,chromeNum);
fitnessPop = fitnessEvaluation(population,f);

% Keep track of the fittest chromosome
fittest = zeros(generationsNum,1);
fittest(1) = max(fitnessPop);
fprintf("Generation 1 \nCurrent Fittest Evaluation: %f \n",fittest(1))

% Crossover decay through generations
decayFun = @(mag,gen) mag*exp(2*(1-gen)/generationsNum);

%% Genetic Algorithm

% Count number of generations
generations = 1;

while generations <= generationsNum
    
    % Keep the previous population
    prevPopulation = population;
    
    % Selection
    population = selectionProcess(population,fitnessPop);
    
    % Crossover - Crossover Parameter = 80%
    population = crossoverProcess(population,0.80);%decayFun(0.80,generations));
    
    % Mutation - Mutation Parameter = 20%
    population = mutationProcess(population,0.20);%decayFun(0.01,generations));
    
    % Elitism - Rate of Chromosomes to be passed directly = 5%
    population = elitismProcess(population,prevPopulation,fitnessPop,0.05);
    
    % Fitness Evaluation
    fitnessPop = fitnessEvaluation(population,f);
    
    % Print Progress Message
    fittest(generations) = max(fitnessPop);
    if mod(generations,generationsNum/10) == 0
        fprintf("Generation %d \nCurrent Fittest Evaluation: %f \n",...
            generations,fittest(generations))
    end
    
    % Count number of generations
    generations = generations + 1;
    
end

%% Results

% Best chromosome of the final population
[~,bestIdx] = max(fitnessPop);
optimalChromosome = population(:,:,bestIdx);
fprintf("\nFitness of the best chromosome: %f \n", fitnessPop(bestIdx) )

% Plot Function Estimation
f_hat = @(u) ObjectiveFuncEstim(u,optimalChromosome);
figure
plot(U,arrayfun(f_hat,U))
title("Estimated Function")
axis tight

% Plot Fittest Chromosome's history
figure
plot(1:generationsNum,fittest)
title('Fittest Chromosome - Fitness Evaluation through Generations')
xlabel('Generations')
ylabel('Fitness Evaluation')

U1 = [-2,0.0012,2];
MSE = sum((arrayfun(f_hat,U1)-arrayfun(f,U1)).^2)/length(U1);
fprintf("Mean Square Error: %f \n",MSE);

toc

%% Functions

%% Initialize Population

function chrom = initPopulation(gaussiansNumber,chromeNumber)

% 3 Genes for each Gaussian
% Magnitude, Center, Std
geneNumber = 3;

% Random Initialize Population of Chromosomes
% Every chromosome includes several genes (gaussians)
% All Parameters initialized in [-4,4]
chrom = rand([gaussiansNumber,geneNumber,chromeNumber])*8-4;

end

%% Fittness Evaulation

% terminate the Genetic Algorithm after a specific event
function fitChrome = fitnessEvaluation(population,f)

% Dataset to be used
U = -2:0.001:2;

fitChrome = zeros(size(population,3),1);

% For every input set of inputs (u1,u2) in the Dataset calculate fitness of the chromosome i
for i = 1:size(population,3)
    for u1 = U
        fitChrome(i) = fitChrome(i) + fitnessCalc(f,[u1],population(:,:,i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fitChrome = fitChrome/length(U);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%% Selection Proccess
% Roulette Wheel Selection
function newPopulation = selectionProcess(population,fitnessPop)

newPopulation = zeros(size(population));

counter = 1;

% Choose Chromosomes to survive
while counter <= size(newPopulation,3)
    
    % Split probabilities into segments
    % Ex: Segment 1 : [0,prob(1)]
    %     Segment 2 : [prob(1),prob(2)] etc.
    prob = length(fitnessPop);
    for i = 1:length(fitnessPop)
        prob(i) = sum(fitnessPop(1:i))/sum(fitnessPop);
    end
    
    % Generate a random number
    randNum = rand;
    
    % Choose which chromosome survives depending
    % on which segment includes the random number
    for index = 1 : length(fitnessPop)
        if prob(index)>=randNum
            break;
        end
    end
    
    % Pass the "lucky" chromosome to the next step
    newPopulation(:,:,counter) = population(:,:,index);
    
    counter = counter + 1;
    
end

end

%% Crossover Proccess
function newPopulation = crossoverProcess(population,crossParam)


offspring = population;

counter = 1;

while counter <= size(offspring,3)
    
    prob = rand;
    
    if prob < crossParam
        
        % Choose 2 parent-chromosomes randomly
        parents = population(:,:,randi([1,size(population,3)],2,1));
        
        % Choose position of gene (gaussian) exchange randomly
        crossPos = randi([1,size(population,1)-1]);
        
        % ie. parent1   = [gene11; gene12; gene13]
        %     parent2   = [gene21; gene22; gene23]
        %     crossPos  = 2
        %     offspring = [gene11; gene12; gene23]
        offspring(:,:,counter) = [parents(1:crossPos,:,1);parents(crossPos+1:end,:,2)];
        
    end
    
    counter = counter + 1;
    
end

% Shuffle the population
idx = randperm(size(population,3));
newPopulation = offspring(:,:,idx);

end


%% Mutation Process
function newPopulation = mutationProcess(population,mutParam)

newPopulation = population;

counter = 1;

while counter <= size(population,3)
    
    prob = rand;
    
    if prob < mutParam
        
        % Random Gene's Position
        randGenePos = randi([1,size(population,2)]);
        
        % New Random Gene ~ Normal Distribution
        newGene = randn(size(population,2),1)*2;
        
        % Random Gene Mutation
        newPopulation(randGenePos,:,counter) = newGene';
    end
    
    counter = counter + 1;
    
end


end

%% Elitism Process
function newPopulation = elitismProcess(population,prevPopulation,fitnessPop,elitParam)

newPopulation = population;
[~,Idx] = sort(fitnessPop,'descend');
for i = 1 : ceil(elitParam*size(population,3))
    newPopulation(:,:,Idx(end+1-i))= prevPopulation(:,:,Idx(i));
end

end


%% Other Functions

% Fitness Function
%[Objective Function: f, Input Vector: u, Parameter's Matrix: chromosome]
function fit = fitnessCalc(f,u,chromosome)

% Choose fittness function formula - normalized
% fitnessFunc = @(x) 1 / (1+abs(x));
fitnessFunc = @(x) 1 / (1+x.^2);

f_hat = ObjectiveFuncEstim(u,chromosome);

% Substruct target function to get the error
error = f_hat - f(u);

fit = fitnessFunc(error);

end

% Objective Function Estimation
% [Input: u , Parameter's Matrix: chromosome]
function func = ObjectiveFuncEstim(u,chromosome)

% Every row of has 3 genes - gaussian Parameters
G = @(gauss) gauss(1)*exp(-(u(1)-gauss(2))^2/(2*gauss(3)^2));

func = 0;
% Add all Gaussians
for i = 1:size(chromosome,1)
    func = func + G(chromosome(i,:));
end

end