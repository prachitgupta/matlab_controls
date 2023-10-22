% doc ctrb
% doc canon
% doc obsv
% doc rank
% doc ss
%doc fprintf

%%Q1
A1 = [-1 0; 0 -2];
B1 = [1 ; sqrt(2)];
C1 = [1 -sqrt(2)/2];
D1 = 0;

is_controlable1 = access_controllability(A1,B1);
is_observable1 = access_observability(A1,C1);
[A_c1,B_c1,C_c1, D_c1,A_o1,B_o1,C_o1, D_o1] = Compute_canonical_forms(A1,B1,C1,D1);

%%%Q2
A2 = [0 1 0; 0 0 1 ; -52 -30 -4];
B2 = [0;0;1];
C2 = [20 1 0];
D2 = 0;

is_controlable2 = access_controllability(A2,B2);
is_observable2 = access_observability(A2,C2);
[A_c2,B_c2,C_c2, D_c2,A_o2,B_o2,C_o2, D_o2] = Compute_canonical_forms(A2,B2,C2,D2);


%%%Q3
A3 = [0 1 0 0; 0 0 1 0; 0 0 0 1; -962 -126 -67 -4];
B3 = [0; 0; 0; 1];
C3 = [300 0 0 0];
D3 = 0;

is_controlable3 = access_controllability(A3,B3);
is_observable3 = access_observability(A3,C3);
[A_c3,B_c3,C_c3, D_c3,A_o3,B_o3,C_o3, D_o3] = Compute_canonical_forms(A3,B3,C3,D3);



%%Q4
A4 = [0 1 0 0; 0 0 1 0; 0 0 0 1; -680 -176 -86 -6];
B4 = [0; 0; 0; 1];
C4 = [100 20 10 0];
D4 = 0;

is_controlable4 = access_controllability(A4,B4);
is_observable4 = access_observability(A4,C4);
[A_c4,B_c4,C_c4, D_c4,A_o4,B_o4,C_o4, D_o4] = Compute_canonical_forms(A4,B4,C4,D4);


%%%functions
function is_controlable = access_controllability(A,B)
%%access system controllability
Con_Matrix = ctrb(A,B);
is_controlable = false;
if(rank(Con_Matrix) == size(A,1))
    is_controlable = true;
    disp("controllable as full rank");
else
    fprintf("Uncontrollable as rank = %f < %f \n",rank(A,1),size(A,1));
end
disp("controllability matrix = ")
disp(Con_Matrix);
end

function is_observable = access_observability(A,C)
%%%access system observability
Obs_Matrix = obsv(A,C);
is_observable = false;
if(rank(Obs_Matrix) == size(A,1))
    is_observable = true;
    disp("observable as full rank");
else
    fprintf("Unobservable as rank = %f < %f \n",rank(A,1),size(A,1));
end
disp("observability matrix = ")
disp(Obs_Matrix);
end

%%%Compute the canonical forms.
function [A_c,B_c,C_c, D_c, A_o,B_o,C_o, D_o] = Compute_canonical_forms(A,B,C,D)
sys = ss(A,B,C,D);
[Obs_canonical_form,T1] = canon(sys,"companion");
A_o = Obs_canonical_form.A;
B_o = Obs_canonical_form.B;
C_o = Obs_canonical_form.C;
D_o = Obs_canonical_form.D;
A_c = A_o.';
B_c = C_o.';
C_c = B_o.';
D_c = D_o;
disp('A matrix of Controller Canonical Form:');
disp(A_c);

disp('B matrix of Controller Canonical Form:');
disp(B_c);

disp('C matrix of Controller Canonical Form:');
disp(C_c);

disp('D matrix of Controller Canonical Form:');
disp(D_c);

disp('A matrix of Observer Canonical Form:');
disp(A_o);

disp('B matrix of Observer Canonical Form:');
disp(B_o);

disp('C matrix of Observer Canonical Form:');
disp(C_o);

disp('D matrix of Observer Canonical Form:');
disp(D_o);
end





