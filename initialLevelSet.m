function Psi = initialLevelSet(circleRadius, sz_i, sz_j)

disp( 'Click for center of initial contour' ); % Contour is initialized manualy..

centers = round(ginput( 1 )); 
c1 = centers(1, :);
Psi = zeros(sz_i, sz_j);
c = sqrt(2) / 2;
s = sqrt(1 - c ^ 2);
for i = 1:sz_i
    for j = 1:sz_j
    ic = (i-c1(2));
    jc = (j-c1(1));
    if( sqrt( ((ic*c+jc*s) / 2.0)^2 + (-jc*s+ic*c)^2) < circleRadius ) 
        Psi(i,j)=1;
    end
    end
end
Psi = -2 * Psi + 1;
%Psi = generateLevelSet(Psi);
end