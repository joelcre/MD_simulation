%----main program for assignment 1.8-----%

clear all
clc
close all
%filename1 = 'energy.txt';
filename2 ='dist_0.txt';
filename3 ='dist_1.txt';
filename4 ='dist_2.txt';
filename5 ='dist_3.txt';
filename6 ='dist_4.txt';
filename7 ='dist_5.txt';
%filename5='distances.txt';
%A1 = importdata(filename1);
%dist=importdata(filename5);
dim1=importdata(filename2);

x=dim1(:,1);
y=dim1(:,2);
z=dim1(:,3);


dim2=importdata(filename3);

x2=dim2(:,1);
y2=dim2(:,2);
z2=dim2(:,3);
dim3=importdata(filename4);

x3=dim3(:,1);
y3=dim3(:,2);
z3=dim3(:,3);

dim4=importdata(filename5);

x4=dim4(:,1);
y4=dim4(:,2);
z4=dim4(:,3);

dim5=importdata(filename6);

x5=dim5(:,1);
y5=dim5(:,2);
z5=dim5(:,3);

dim6=importdata(filename7);

x6=dim6(:,1);
y6=dim6(:,2);
z6=dim6(:,3);

%X=dist(:,1);
%Y=dist(:,2);
%Z=dist(:,3);






%pot_en= A1(:,2);
%kin_en=A1(:,3);
%tot_en=A1(:,1);

% cut off the first values because the simulation has not reached
% equilibrium
%pot_en_eq=pot_en(401:1:length(pot_en));
%kin_en_eq=kin_en(401:1:length(kin_en));
%tot_en_eq=tot_en(401:1:length(tot_en));

%%figure(1)
%%plot(pot_en,'b')
%%hold on
%%plot(kin_en,'r') 
%%hold on
%%plot(tot_en,'g')
%
%xlabel('time')
%ylabel('energy')

figure(1)
%mesh(dim,'.')
plot3(x,y,z,'.b')
hold on
plot3(x2,y2,z2,'.r')
hold on
plot3(x3,y3,z3,'.g')
hold on
plot3(x4,y4,z4,'.y')
hold on
plot3(x5,y5,z5,'.c')
hold on
plot3(x6,y6,z6,'.k')

xlabel('x')
ylabel('y')
zlabel('z')
