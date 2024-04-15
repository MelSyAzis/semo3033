clear all;
clc;

L = 1; % length
N = 5; % no of nodes
dx = L / N; % CV length
x = [dx/2 : dx : L]; # node coordinates
TA = 0;
TB = 10000;

A = [
300 -100 0 0 0;
-100 200 -100 0 0;
0 -100 200 -100 0;
0 0 -100 200 -100;
0 0 0 -100 300
]

b = [
200 * TA;
0;
0;
0;
200 * TB;
]

T = A \ b;

plot(x, T, '*-')

xlabel('x')

ylabel('T')