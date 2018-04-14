%----main program for assignment 1.8-----%

clear all
clc
close all
filename1 = 'energy.txt';
A1 = importdata(filename1);

pot_en= A1(:,2);
kin_en=A1(:,3);
tot_en=A1(:,1);

% cut off the first values because the simulation has not reached
% equilibrium
%pot_en_eq=pot_en(401:1:length(pot_en));
%kin_en_eq=kin_en(401:1:length(kin_en));
%tot_en_eq=tot_en(401:1:length(tot_en));


plot(pot_en,'b')
hold on
plot(kin_en,'r') 
hold on
plot(tot_en,'g')

xlabel('time')
ylabel('energy')




