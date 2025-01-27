clear all
clc

engineSpan = 21.1; % distance from centreline to outer engine (m)
engineToTail = 36; % distance from engine to vstab (m)
T = 321600; % max. thrust (N)
sideF = T * engineSpan / engineToTail; % required vstab sideforce (N)