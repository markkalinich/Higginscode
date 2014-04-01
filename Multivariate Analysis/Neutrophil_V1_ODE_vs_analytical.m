%This code was written by Mark Kalinich on 20140305 to use ODE45 to
%numerically integrate a system of equations describing neutrophil
%population kinetics.  

%Differential Equation Model: 

%d(New)/dt = R_production,new - R_death,old - R_maturation
%d(Old)/dt = R_production,old - R_death,old + R_maturation

%for Steady state: 0=R_production,old - R_death,old + R_maturation
%0 = (Cell Production rate)*(Old Fraction)-(Total Death)*(Old/(Old+New)
%for only old cells:
%0 = Cell Production Rate - Total Death Rate

%For Only Young Cells: 

%Model Assumptions: 
%1) Constant rate of total cell production from bone
%2) constant ratio of old to new neutrophils in the bone
%3) random selection of neutrophils from bone
%4) Constant fraction of young neutrophils mature into old neutrophils
%5) Constant rate of total cell death (Death_old = Constant*(# dead cells)
%6) Cell death is not dependent on cell age
%7) R_production is a binomial random variable that uses the specified
%fraction of old vs. new cells. (R_production_total= constant = R_new+R+old) 
%8) Person weighs 70 kg. 

clc;close all;clear all;

%set integration boundaries for weight of catalyst
T0=0; % minutes
TF=2000; %minutes, ~4 days

%set initial conditions, making sure units are compatible
Weight                = 70; %person's weight, in kg 
Young_fraction        = 0.1;
Neutrophil_blood      = 65*10^7*Weight; %neutrophils initially in blood pool for a 70 kg person, "Neutrophil kinetics in health and disease" Summers 2010 
Band_Fraction         = .015; % For an uninfected person, bands range from 0-3% in blood http://www.nlm.nih.gov/medlineplus/ency/article/003657.htm
Young_initial         = Neutrophil_blood*Band_Fraction; % initial number of bands in the pool 
Old_initial           = Neutrophil_blood*(1-Band_Fraction);% initial number of mature cells in pool
Cells_Entering_System = 1.7*10^9*Weight;  %this is the number of cells entering from bone marrow per day "Neutrophil kinetics in health and disease" Summers 2010
Mature_probability    = 1; %fraction of band cells in blood that will mature after 1 day
Death_Leaving_System  = 1.7*10^9*Weight; %assume steady state
Mean_residence_time   = Neutrophil_blood/Cells_Entering_System;
%Add in conversions to get proper rates.
Time_Convert = 1440; %minutes/day
Cell_Convert = 1*10^6; %convert cells to millions of cells
Cell_Enter   = round(Cells_Entering_System/(Time_Convert*Cell_Convert)); 
Mature_Prob  = Mature_probability/Time_Convert; %fraction of band cells that  matures into an old cell in 1 minute.
Death        = round(Death_Leaving_System/(Time_Convert*Cell_Convert));% number of cells to die each cycle; assume steady state

%pack initial conditions
N0=[Young_initial/Cell_Convert Old_initial/Cell_Convert]; %10^6 cells

%pack constant parameters
p.Y_enter   = Cell_Enter*Young_fraction; %cells/minute
p.Old_enter = Cell_Enter*(1-Young_fraction); %K
p.Mature    = Mature_Prob; %mature probability in minutes
p.Death     = Death; %1/Pa


%execute the ode solver
[T,N]=ode45(@(T,N)odeneutrophil(T,N,p),[T0 TF],N0);

%Compare to Analytical Solution
k            = Death./(Neutrophil_blood./Cell_Convert)+Mature_Prob; %exp(-kt)f or model
Young_SS     = Cell_Enter*Young_fraction/k;


%Calculate Young and Old Neutrophil Counts
Time = T0:TF; %
Young = (Young_initial/Cell_Convert).*exp(-k.*Time)+(Cell_Enter*Young_fraction/k)*(1-exp(-k.*Time));
Old   = Neutrophil_blood/Cell_Convert - Young;



fig = figure;
plot(Time,Young,'gd');
hold on
plot(T,N(:,1),'ro');
hold on
plot(Time,Young_SS,'b')

legend('Analytical','Numerical','Steady-State')
xlabel ('Time (minutes)')
ylabel ('Young Neutrophils (10^6 cells)')
title ('Analytical vs.Numerical Young Neutrophil Count')


% fig = figure;
% plot(T,N(:,1),'r');
% hold on
% plot(T, N(:,2),'g');
% hold on
% plot(T, N(:,1)+N(:,2), 'b');
% xlabel ('time (minutes)')
% ylabel ('Neutrophils')
% %title (strcat('Neutrophils vs. Time, Const. Tot. Death Rate, Young Neutrophil Mature Probability (per day) =',num2str(Mature_probability)))
% legend ('Young Neutrophils','Old Neutrophils','Total Neutrophils')
