%Initial 1 Compartment Model of Neutrophil Maturation and Death.  The model 
%starts by adding a set number of cells to a blood compartment, allows some
%young neutrophils to mature to old, and then removes a set number of cells 
%from the population (death).  This code was written by Mark Kalinich on 20140225. 

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


clc; clear;
%First, define my constants
Weight                = 70; %person's weight, in kg
Neutrophil_blood      = 65*10^7*Weight; %neutrophils initially in blood pool for a 70 kg person, "Neutrophil kinetics in health and disease" Summers 2010 
Band_Fraction         = 0.015; % For an uninfected person, bands range from 0-3% in blood http://www.nlm.nih.gov/medlineplus/ency/article/003657.htm
Young_initial         = Neutrophil_blood*Band_Fraction; % initial number of bands in the pool 
Old_initial           = Neutrophil_blood*(1-Band_Fraction);% initial number of mature cells in pool
Cells_Entering_System = 1.7*10^9*Weight;  %this is the number of cells entering from bone marrow per day "Neutrophil kinetics in health and disease" Summers 2010 
Young_fraction        = .1; % select Bone neutrophil fraction I'm exploring

Mature_probability    = .2; %fraction of band cells in blood that will mature after 1 day
Death_Leaving_System  = 1.7*10^9*Weight; %assume steady state
Band_Mature           = .1; %fraction of bands that will mature in 1 day in the blood
Mean_residence_time   = Neutrophil_blood/Cells_Entering_System;
                 
%Add in conversions to get proper rates.
Time_Convert = 1440; %minutes/day
Cell_Convert = 1*10^6; %convert cells to millions of cells
Cell_Enter   = round(Cells_Entering_System/(Time_Convert*Cell_Convert)); 
Mature_Prob  = Mature_probability/Time_Convert; %fraction of band cells that  matures into an old cell in 1 minute.
Death        = round(Death_Leaving_System/(Time_Convert*Cell_Convert));% number of cells to die each cycle; assume steady state
k            = Death./(Neutrophil_blood./Cell_Convert+Mature_Prob); %exp(-kt)f or model


%Calculate Young and Old Neutrophil Counts
Time = 1:1440; %make calculation for 1 day
Young = (Young_initial/Cell_Convert).*exp(-k.*Time)+(Cell_Enter*Young_fraction/k)*(1-exp(-k.*Time));
Old   = Neutrophil_blood/Cell_Convert - Young;

% fig = figure;
% plot(time,Neutrophil_Count(1:length(time),2),'r');
% hold on
% plot(time,Neutrophil_Count(1:length(time),1),'g');
% xlabel ('time (minutes)')
% ylabel ('Neutrophils')
% title (strcat('Neutrophils vs. Time, Const. Tot. Death Rate, Young Neutrophil Mature Probability (per day) =',num2str(Mature_probability)))
% legend ('Mature Neutrophils','Band Neutrophils')

%saveas(fig,strcat('Neutrophil_Kinetics_Young_Neutrophil_Mature_Prob=_',num2str(Mature_probability)),'jpeg')



%end



