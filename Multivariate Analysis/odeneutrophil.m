function dNdT = odeneutrophil(T,N,p)
%function calculates derivatives of species in a PBR model of a catalytic
%converter



%unpack parameters and variables
Young = N(1);% number of young neutrophils
Old   = N(2);% initial number of old neutrophils

%unpack parameters
Y_enter     = p.Y_enter; %rate of young neutrophil entry
Old_enter   = p.Old_enter; %rate of old neutrophil entry
Mature_rate = p.Mature;   %mature rate of young neutrophil to old 
Death       = p.Death;    %death rate 



%calculate derivative of NO using provided rate law
dNdT(1) = Y_enter   - Death*(Young/(Young+Old)) - Mature_rate*Young;
dNdT(2) = Old_enter - Death*(Old/(Young+Old))   + Mature_rate*Young;



%ode solver requires a column vector
dNdT=dNdT';
end