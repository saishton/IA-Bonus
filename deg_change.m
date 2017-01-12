function changes = deg_change(DinT)

M1 = DinT(:,1:end-1);
M2 = DinT(:,2:end);

changes = M2-M1;

end