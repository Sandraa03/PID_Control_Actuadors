//Limpiamos consola
clear
clc

//Declaración de variable

thau=0.001;
gammA=20;
//[theta]=linspace(1,10,10);
//theta=theta';

//Declaración de la kp

kp=10;

//Declaración del vector Kd con valores del 1 al 10 de 1 en 1 para el control

[kd]=linspace(1,10,10);

//Ecuación de control

(theta.^3)-(2*(theta.^2))+(1+(gammA*kd*thau.^2)+(gammA*kd*thau))*(theta.^2)-((gammA*kd*thau)*theta)==0; //theta(s) son las frecuencias naturales del sistema


//Plot de theta en función de Kd

scf(0);
plot2d(theta,kd);
