% Levent Batakci
% Super basic testing for my own sanity

a = 1:5;
b = 1:2:11;

T1 = toeplitzmat(a,b);
T2 = T1 + T1;
T3 = 2.* T1
T6 = -T1
T4 = T3 ./ 2
T5 = T3 ./ T3


T1 >= 5
T1 <= 5
T1 > 5
T1 < 5
T1 ~= 5
T1 == 5

T1 >= 5 & T1 <= 5
T1 >= 5 | T1 <= 5
~(T1 >= 5)

disp("last test")
(T1 > 1) < (1>0)
(T1 > 1) & (1>0)

T1'
T1

A = T1;
B = toeplitzmat([1 flip(A.tr(3:end))], 1:10);
[A B]