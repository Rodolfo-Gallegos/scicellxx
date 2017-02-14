%% Herramienta de Evaluacion de resultados navegaci�n Inercial
%  Desarrollador X(MAJ)

clear
close all
clc
M= load('navigation_data_from_GPS.dat');
clf
%% CARGAR ARCHIVO DE POSICIONES
t   =  M(:,1);          % estampa de tiempo (seg)

lon = -M(:,2);          % longitud(� decimales) 
                        % correcci�n de signo para (East+)

lat =  M(:,3);          % latitud(� decimales) 

distancia= M(:,7);      % distancia relativa entre puntos (m)
var_rumbo= -M(:,8);       % rumbo del veh�culo (�) 
                          % correcci�n de signo 
%% CALCULO

rumbo_inicial=122;  %Rumbo inicial(�) desde el cual comienza la nevegaci�n 
[new_lat,new_lon,Error]=calc_nav(lat,lon,distancia,rumbo_inicial,var_rumbo);

%% RESULTADOS
plot(lon,lat,'r')
hold on
plot(new_lon,new_lat,'.g','MarkerSize',5)
plot_google_map2('Refresh','1','maptype','terrain','AutoAxis','1')

figure
plot(t,Error)
grid on
xlabel('tiempo(seg)')
ylabel('Error(m)')





