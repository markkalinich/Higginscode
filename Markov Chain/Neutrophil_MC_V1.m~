%This code was written to solve a simplifed homogenous discrete-time
%Markov-chain model of neutrophil maturation.  It was written by Mark
%Kalinich on 20140409.  

%Current assumptions
%All Cells in a compartment mature from 1 compartment to the next every
%time step
%BM cell age is equally distributed across ages. 

clc;close all;clear all;
%Step 1: Define my time intervals
t0 = 1; %initial time step
tf = 5; %right now, we're doing things in minutes
Weight                = 70; %person's weight, in kg 
Neutrophil_blood      = 65*10^7*Weight; %neutrophils initially in blood pool for a 70 kg person, "Neutrophil kinetics in health and disease" Summers 2010 
Band_Fraction         = .015; % For an uninfected person, bands range from 0-3% in blood http://www.nlm.nih.gov/medlineplus/ency/article/003657.htm
Young_initial         = Neutrophil_blood*Band_Fraction; % initial number of bands in the pool 
Old_initial           = Neutrophil_blood*(1-Band_Fraction);% initial number of mature cells in pool
Cells_Entering_System = 1.7*10^9*Weight;  %this is the number of cells entering from bone marrow per day "Neutrophil kinetics in health and disease" Summers 2010


%Convert to more easily-usable numbers
N_blood = Neutrophil_blood/10^6; % all counts are now "million cells")
Cell_Enter_Total = Cells_Entering_System/(60*24*10^6); %now we have total number of cells (millions) entering in per day

%start up vectors/matrices
n0 = zeros(tf+1,1); %initialize the n0 vector (intitial value of all ages)assume neutrophils are in first compartment to start
n0(1) = 1; %need to multiply by the BM cell injection rate
n0(2) = N_blood; %initially assume that we start with all cells in youngest compartment
BM = ones(1,tf); %initialize BM injection vector
BM = BM*(Cell_Enter_Total/tf); %assumes that the rate of cell ejection from the BM is equivalent accross ages  
Mature = ones(tf,1); %all cells from a given compartment mature to the next compartment 
Mature(tf) = 0; %cells in last compartment do not mature
Death  = zeros(tf,1); %assume all cells except for the last compartment do not die
Death(tf) = 1; %last compartment dies
alpha = 1-(Mature+Death); %term for 

A = diag(alpha,1); %initialize parameter matrix
A = A([1:tf],[1:tf+1]);
A(:,1) = BM';
for i = 2:tf
A(i,i) = Mature(i-1);
end 
%now we have our A matrix determined for the system.  We need to iterate
%through time to get the number of cells in each compartment as a function
%of time. 
%initialize solution matrix
S = zeros(tf,tf);
S(1,1) = N_blood;
for i = 1:tf
    for j = 1:tf
    n(j,i)= A*(1
 




