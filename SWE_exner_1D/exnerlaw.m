function [ Qu ] = exnerlaw( u,hs,t )
% Exner law (grass law).
% alfa=1/.65;
% Qu=alfa.*sign(u).*abs(u).^1.5;

Af=1.35;
Qu=Af*hs.*u.^3;
end