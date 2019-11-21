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
f = @(u) sin(u(1)+u(2))*sin(u(1).^2);

%% Plot Function to be etsimated

% Count Number of Plots
plotNum = 0;

plotNum = plotFunction(plotNum,f);

figure(plotNum-1)
title('3D Plot - $$ f(u_1,u_2) = sin(u_1+u_2) \cdot sin(u_1^2) $$','Interpreter','Latex')
xlabel('$$ u_1 $$','Interpreter','Latex')
ylabel('$$ u_2 $$','Interpreter','Latex')

figure(plotNum)
title('2D Plot - $$ f(u_1,u_2) = sin(u_1+u_2) \cdot sin(u_1^2) $$','Interpreter','Latex')
xlabel('$$ u_1 $$','Interpreter','Latex')
ylabel('$$ u_2 $$','Interpreter','Latex')

pause(0.01);

%% Initialization

% Number of Chromosomes in every Population
chromeNum = 50;

% Number of Gaussians in every Chromosome
gaussiansNum = 15;

% Number of generations until termination
generationsNum = 100;

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
    
    % Selection
    population = selectionProcess(population,fitnessPop);
    
    % Crossover - Initial Crossover Parameter = 50%
    population = crossoverProcess(population,0.50);%decayFun(0.50,generations));
    
    % Mutation - Initial Mutation Parameter = 30%
    population = mutationProcess(population,0.30);%decayFun(0.30,generations));
    
    % Fitness Evaluation
    fitnessPop = fitnessEvaluation(population,f);
    
    % Fittest Chromosome's Evaluation
    fittest(generations) = max(fitnessPop);
    
    if mod(generations,generationsNum/10) == 0
        fprintf("Generation %d \nCurrent Fittest Evaluation: %f \n",...
            generations,fittest(generations))
    end
    
    % Count number of generations
    generations = generations + 1;
    
end

%% Results - Evaluation


% Best chromosome of the final population
[~,bestIdx] = max(fitnessPop);
optimalChromesome = population(:,:,bestIdx);
fprintf("\nFitness of the best chromosome: %f \n", fitnessPop(bestIdx) )

f_hat = @(u) ObjectiveFuncEstim(u,optimalChromesome);

% Mean Square Error
MSE = 0;
U_test = -2:0.015:2;

for u1 = U_test
    for u2 = U_test
        MSE = MSE + (f_hat([u1,u2])- f([u1,u2])).^2/length(U_test)^2;
    end
end
fprintf("\nMean Square Error: %f \n",MSE)

% Plot Function Estimation
plotNum = plotFunction(plotNum,f_hat);

figure(plotNum-1)
title('3D Plot - Estimated $$ \hat{f}(u_1,u_2)  $$','Interpreter','Latex')
xlabel('$$ u_1 $$','Interpreter','Latex')
ylabel('$$ u_2 $$','Interpreter','Latex')

figure(plotNum)
title('3D Plot - Estimated $$ \hat{f}(u_1,u_2)  $$','Interpreter','Latex')
xlabel('$$ u_1 $$','Interpreter','Latex')
ylabel('$$ u_2 $$','Interpreter','Latex')

plotNum = plotNum + 1;

% Plot Fittest Chromosome's history
figure(plotNum)
plot(1:generationsNum,fittest)
title('Fittest Chromosome - Fitness Evaluation through Generations')
xlabel('Generations')
ylabel('Fitness Evaluation')

toc

%% Functions

%% Initialize Population

function chrom = initPopulation(gaussiansNumber,chromeNumber)

% 5 Genes for each Gaussian
% Magnitude, Center1, Center2, Std1, Std2
geneNumber = 5;

% Random Initialize Population of Chromosomes
% Every chromosome includes several genes (gaussians)
% All Parameters initialized in [-4,4]
chrom = rand(gaussiansNumber,geneNumber,chromeNumber)*8-4;

end

%% Fittness Evaulation

% terminate the Genetic Algorithm after a specific event
function fitChrome = fitnessEvaluation(population,f)

% Dataset to be used
U = -2:0.05:2;

fitChrome = zeros(size(population,3),1);

% For every input set of inputs (u1,u2) in the Dataset calculate fitness of the chromosome i
for i = 1:size(population,3)
    
    for u1 = U
        for u2 = U
            fitChrome(i) = fitChrome(i) + fitnessCalc(f,[u1;u2],population(:,:,i));
        end
    end
    
end

fitChrome = fitChrome/length(U)^2;

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
    
    % Pass the "lucky" chromosome to the next generation
    newPopulation(:,:,counter) = population(:,:,index);
    
    counter = counter + 1;
    
end

end

%% Crossover Proccess
function newPopulation = crossoverProcess(population,crossParam)


offspring = zeros(size(population,1),size(population,2),ceil(size(population,3)*crossParam));

counter = 1;

while counter <= size(offspring,3)
    
    % Choose 2 parent-chromosomes randomly
    parents = population(:,:,randi([1,size(population,3)],2,1));
    
    % Choose position of gene (gaussian) exchange randomly
    crossPos = randi([1,size(population,1)-1]);
    
    % ie. parent1   = [gene11; gene12; gene13]
    %     parent2   = [gene21; gene22; gene23]
    %     crossPos  = 2
    %     offspring = [gene11; gene12; gene23]
    offspring(:,:,counter) = [parents(1:crossPos,:,1);parents(crossPos+1:end,:,2)];
    
    counter = counter + 1;
    
end

% Join offsprings with part of the elder population
newPopulation = population;
newPopulation(:,:,1:size(offspring,3)) = offspring;

% shuffle the population
idx = randperm(size(newPopulation,3));
newPopulation = newPopulation(:,:,idx);

end


%% Mutation Proccess
function newPopulation = mutationProcess(population,mutParam)

newPopulation = population;

counter = 1;

% Random number of mutations take place
mutChance = 2*mutParam*rand;

while counter <= round(size(newPopulation,3)*mutChance)
    
    % Chromosome Mutation
%     newPopulation(:,:,randi([1,size(population,3)])) = ...
%         randn(size(population,1),size(population,2),1)*8-4;

    % Gene Mutation of a random Chromosome
    newPopulation(:,randi([1,size(population,2)]),randi([1,size(population,3)])) = ...
        randn(size(population,1),1,1)*8-4;

    counter = counter + 1;
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

% Every row of has 5 genes - gaussian Parameters
G = @(gauss) gauss(1)*exp(-(u(1)-gauss(2))^2/(2*gauss(4)^2)-(u(2)-gauss(3))^2/(2*gauss(5)^2));

func = 0;
% Add all Gaussians
for i = 1:size(chromosome,1)
    func = func + G(chromosome(i,:));
end

end

% Plot the function f (3D + 2D)
function plotNum = plotFunction(plotNum,f)

x = linspace(-2, 2, 100);
y = linspace(-2, 2, 100);

[X,Y] = meshgrid(x,y);

func = [];
for i = 1:length(x)
    for j = 1:length(y)
        func(i,j) = f([X(i,j);Y(i,j)]);
    end
end

figure(plotNum+1)
surf(X,Y,func)
view(160,40)
colorbar

figure(plotNum + 2)
contour(X,Y,func,20)
colorbar

plotNum = plotNum + 2;

end

% Function to automatically save plots in high resolution
function savePlot(name)

% Resize current figure to fullscreen for higher resolution image
set(gcf, 'Position', get(0, 'Screensize'));

% Save current figure with the specified name
saveas(gcf, join([name,'.jpg']));

% Resize current figure back to normal
set(gcf,'position',get(0,'defaultfigureposition'));

end