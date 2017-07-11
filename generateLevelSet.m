% Psi_ : inside -1, outside 1
function Psi = generateLevelSet(Psi_)

Psi = double((Psi_ > 0).*(bwdist(Psi_ < 0) - 0.5) - (Psi_ < 0).*(bwdist(Psi_ > 0) - 0.5));

end