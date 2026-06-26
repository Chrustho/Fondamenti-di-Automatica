%% Esercizio C - Catena di Markov a tempo discreto
%  Fondamenti di Automatica - Traccia 28 - Bruni Christian, matr. 240008
%
%  Sorgente equivalente al Live Script ParametriMarkov.mlx.
%  Definisce i parametri (matrice di transizione e stato iniziale) che
%  alimentano il modello Simulink Simulink2a.slx, dove la ricorsione
%  x(k+1) = A x(k) viene iterata fino alla convergenza allo stato
%  stazionario della catena.

clear; close all;

% Matrice di transizione A. Le colonne sono distribuzioni di probabilità
% (somma = 1): A(i,j) è la probabilità di passare dallo stato j allo
% stato i. Essendo stocastica per colonne, A ammette 1 come autovalore e
% un autovettore associato che è la distribuzione stazionaria.
A = [1/5  5/19  3/11;
     1/3  9/19  7/22;
     7/15 5/19  9/22]

%% Scelta di uno stato iniziale x0 random
% x0 è una distribuzione di probabilità iniziale generata in modo
% pseudo-casuale e poi normalizzata (somma = 1). La scelta è arbitraria:
% per la convergenza conta solo che sia un vettore stocastico.
x0 = rand(3,1)
x0 = x0/sum(x0)
s  = sum(x0)            % verifica: deve valere 1
