% Konstantinos Letros 8851
% Optimization Techniques
% The Project - Parameters Estimation
% Parameters Estimation using Genetic Algorithm

%% Clean the screen

clc
clear
close all;
format long;

%% Plot Function to be etsimated

% Count Number of Plots
plotNum = 0;

plotNum = plotFunction(plotNum);

%% Initialization

% Number of Chromosomes in every Population
chromeNum = 4;

% Number of Gaussians in every Chromosome
gaussiansNum = 15;

% Initialize Population of Chromosomes
% population(gaussian_i,parameter_j,chromosome_k)
population = initPopulation(gaussiansNum,chromeNum);
fitnessPop = fitnessEvaulation(population)

%% Genetic Algorithm

% Termination Criterion
terminate = false;

% Main Loop
while terminate == false
    
    selectionProccess()
    crossoverProccess()
    mutationProccess()
    fitnessEvaluation()
    
end

%% Initialize Population

function chrom = initPopulation(gaussiansNumber,chromeNumber)

% 5 Genes for each Gaussian
% Magnitude, Center1, Center2, Std1, Std2
geneNumber = 5;

% Random Initial Population of Chromosomes
% Every chromosome includes several gaussians
chrom = rand(gaussiansNumber,geneNumber,chromeNumber);

end

%% Fittness Evaulation

% terminate the Genetic Algorithm after a specific event
function fitChrome = fitnessEvaulation(population)

% Set of Data to be used
U = -2:0.1:2;

fitChrome = zeros(size(population,3),1);
for i = 1:size(population,3)
    for u1 = U
        for u2 = U
            fitChrome(i) = fitChrome(i) + fitnessCalc([u1;u2],population(:,:,i));
        end
    end
end

fitChrome = fitChrome./length(U)^2;

end

%% Selection Proccess

function selectionProcess()



end

%% Crossover Proccess

function crossoverProcess()



end


%% Mutation Proccess

function mutationProcess()



end

%% Other Functions

% Function to be estimated (Intput Vector: u)
function func = f(u)
func = sin(u(1)+u(2))*sin(u(1)).^2;
end

% Fitness Function (Input Vector: u, Parameter's Matrix: chromosome)
function fit = fitnessCalc(u,chromosome)
    
    % Choose fittness function formula - normalized
    fitnessFunc = @(x) 1 / (1+abs(x));

    % Every row of has 5 genes - gaussian Parameters
    G = @(gauss) gauss(1)*exp(-(u(1)-gauss(2))^2/(2*gauss(4)^2)-(u(2)-gauss(3))^2/(2*gauss(5)^2));
    
    f_hat = 0;
    % Add all Gaussians
    for i = 1:size(chromosome,1)
     f_hat = f_hat + G(chromosome(i,:));
    end
    
    % Substruct target function to get the error 
    error = f_hat - f(u);
    
    fit = fitnessFunc(error);
    
end

% Plot the function f
function plotNum = plotFunction(plotNum)

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
view(-10,25)
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