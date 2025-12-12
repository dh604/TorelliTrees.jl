from admcycles import *

A2A2_k1_integral = 0;
A2A2_T1_term = 0;

A2A2_T1_term += 3 * hodge_integral(2, [2], [2], [], n=1) * hodge_integral(2, [2], [1], [1], n=1);
print('Computed term 1 / 16 of tree 1 / 9');
A2A2_T1_term += 3 * hodge_integral(2, [2], [1], [1], n=1) * hodge_integral(2, [2], [2], [], n=1);
print('Computed term 2 / 16 of tree 1 / 9');
A2A2_T1_term += -4 * hodge_integral(2, [2, 1], [1], [], n=1) * hodge_integral(2, [2], [1], [1], n=1);
print('Computed term 3 / 16 of tree 1 / 9');
A2A2_T1_term += -2 * hodge_integral(2, [2, 1], [0], [1], n=1) * hodge_integral(2, [2], [2], [], n=1);
print('Computed term 4 / 16 of tree 1 / 9');
A2A2_T1_term += 1 * hodge_integral(2, [2, 1, 1], [0], [], n=1) * hodge_integral(2, [2], [1], [1], n=1);
print('Computed term 5 / 16 of tree 1 / 9');
A2A2_T1_term += 2 * hodge_integral(2, [2, 2], [0], [], n=1) * hodge_integral(2, [2], [1], [1], n=1);
print('Computed term 6 / 16 of tree 1 / 9');
A2A2_T1_term += -2 * hodge_integral(2, [2], [2], [], n=1) * hodge_integral(2, [2, 1], [0], [1], n=1);
print('Computed term 7 / 16 of tree 1 / 9');
A2A2_T1_term += -4 * hodge_integral(2, [2], [1], [1], n=1) * hodge_integral(2, [2, 1], [1], [], n=1);
print('Computed term 8 / 16 of tree 1 / 9');
A2A2_T1_term += 3 * hodge_integral(2, [2, 1], [1], [], n=1) * hodge_integral(2, [2, 1], [0], [1], n=1);
print('Computed term 9 / 16 of tree 1 / 9');
A2A2_T1_term += 3 * hodge_integral(2, [2, 1], [0], [1], n=1) * hodge_integral(2, [2, 1], [1], [], n=1);
print('Computed term 10 / 16 of tree 1 / 9');
A2A2_T1_term += -1 * hodge_integral(2, [2, 1, 1], [0], [], n=1) * hodge_integral(2, [2, 1], [0], [1], n=1);
print('Computed term 11 / 16 of tree 1 / 9');
A2A2_T1_term += -2 * hodge_integral(2, [2, 2], [0], [], n=1) * hodge_integral(2, [2, 1], [0], [1], n=1);
print('Computed term 12 / 16 of tree 1 / 9');
A2A2_T1_term += 1 * hodge_integral(2, [2], [1], [1], n=1) * hodge_integral(2, [2, 1, 1], [0], [], n=1);
print('Computed term 13 / 16 of tree 1 / 9');
A2A2_T1_term += -1 * hodge_integral(2, [2, 1], [0], [1], n=1) * hodge_integral(2, [2, 1, 1], [0], [], n=1);
print('Computed term 14 / 16 of tree 1 / 9');
A2A2_T1_term += 2 * hodge_integral(2, [2], [1], [1], n=1) * hodge_integral(2, [2, 2], [0], [], n=1);
print('Computed term 15 / 16 of tree 1 / 9');
A2A2_T1_term += -2 * hodge_integral(2, [2, 1], [0], [1], n=1) * hodge_integral(2, [2, 2], [0], [], n=1);
print('Computed term 16 / 16 of tree 1 / 9');
A2A2_k1_integral += 1/1 * A2A2_T1_term;

A2A2_T2_term = 0;

A2A2_T2_term += 1 * hodge_integral(2, [2], [2, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 1 / 7 of tree 2 / 9');
A2A2_T2_term += 1 * hodge_integral(2, [2], [1, 1], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 2 / 7 of tree 2 / 9');
A2A2_T2_term += 1 * hodge_integral(2, [2], [0, 2], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 3 / 7 of tree 2 / 9');
A2A2_T2_term += -2 * hodge_integral(2, [2, 1], [1, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 4 / 7 of tree 2 / 9');
A2A2_T2_term += -2 * hodge_integral(2, [2, 1], [0, 1], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 5 / 7 of tree 2 / 9');
A2A2_T2_term += 1 * hodge_integral(2, [2, 1, 1], [0, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 6 / 7 of tree 2 / 9');
A2A2_T2_term += 2 * hodge_integral(2, [2, 2], [0, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 7 / 7 of tree 2 / 9');
A2A2_k1_integral += 1/2 * A2A2_T2_term;

A2A2_T3_term = 0;

A2A2_T3_term += 1 * hodge_integral(2, [2], [2, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 1 / 7 of tree 3 / 9');
A2A2_T3_term += 1 * hodge_integral(2, [2], [1, 1], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 2 / 7 of tree 3 / 9');
A2A2_T3_term += 1 * hodge_integral(2, [2], [0, 2], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 3 / 7 of tree 3 / 9');
A2A2_T3_term += -2 * hodge_integral(2, [2, 1], [1, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 4 / 7 of tree 3 / 9');
A2A2_T3_term += -2 * hodge_integral(2, [2, 1], [0, 1], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 5 / 7 of tree 3 / 9');
A2A2_T3_term += 1 * hodge_integral(2, [2, 1, 1], [0, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 6 / 7 of tree 3 / 9');
A2A2_T3_term += 2 * hodge_integral(2, [2, 2], [0, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 7 / 7 of tree 3 / 9');
A2A2_k1_integral += 1/2 * A2A2_T3_term;

A2A2_T4_term = 0;

A2A2_T4_term += 1 * hodge_integral(1, [1], [1, 0], [], n=2) * hodge_integral(1, [1], [0, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 1 / 6 of tree 4 / 9');
A2A2_T4_term += 1 * hodge_integral(1, [1], [0, 1], [], n=2) * hodge_integral(1, [1], [0, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 2 / 6 of tree 4 / 9');
A2A2_T4_term += 1 * hodge_integral(1, [1], [0, 0], [1], n=2) * hodge_integral(1, [1], [1, 0], [], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 3 / 6 of tree 4 / 9');
A2A2_T4_term += 1 * hodge_integral(1, [1], [0, 0], [1], n=2) * hodge_integral(1, [1], [0, 1], [], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 4 / 6 of tree 4 / 9');
A2A2_T4_term += -2 * hodge_integral(1, [1, 1], [0, 0], [], n=2) * hodge_integral(1, [1], [0, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 5 / 6 of tree 4 / 9');
A2A2_T4_term += -2 * hodge_integral(1, [1], [0, 0], [1], n=2) * hodge_integral(1, [1, 1], [0, 0], [], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 6 / 6 of tree 4 / 9');
A2A2_k1_integral += 1/1 * A2A2_T4_term;

A2A2_T5_term = 0;

A2A2_T5_term += -6 * hodge_integral(0, [0], [0, 0, 0], [], n=3) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(2, [2], [1], [1], n=1);
print('Computed term 1 / 2 of tree 5 / 9');
A2A2_T5_term += 6 * hodge_integral(0, [0], [0, 0, 0], [], n=3) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(2, [2, 1], [0], [1], n=1);
print('Computed term 2 / 2 of tree 5 / 9');
A2A2_k1_integral += 1/2 * A2A2_T5_term;

A2A2_T6_term = 0;

A2A2_T6_term += -6 * hodge_integral(0, [0], [0, 0, 0], [], n=3) * hodge_integral(2, [2], [1], [1], n=1) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 1 / 2 of tree 6 / 9');
A2A2_T6_term += 6 * hodge_integral(0, [0], [0, 0, 0], [], n=3) * hodge_integral(2, [2, 1], [0], [1], n=1) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 2 / 2 of tree 6 / 9');
A2A2_k1_integral += 1/2 * A2A2_T6_term;

A2A2_T7_term = 0;

A2A2_T7_term += -3 * hodge_integral(0, [0], [0, 0, 0], [], n=3) * hodge_integral(1, [1], [0, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 1 / 1 of tree 7 / 9');
A2A2_k1_integral += 1/1 * A2A2_T7_term;

A2A2_T8_term = 0;

A2A2_T8_term += -3 * hodge_integral(0, [0], [0, 0, 0], [], n=3) * hodge_integral(1, [1], [0, 0], [1], n=2) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 1 / 1 of tree 8 / 9');
A2A2_k1_integral += 1/1 * A2A2_T8_term;

A2A2_T9_term = 0;

A2A2_T9_term += -6 * hodge_integral(0, [0], [0, 0, 0, 0], [1], n=4) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1) * hodge_integral(1, [1], [0], [], n=1);
print('Computed term 1 / 1 of tree 9 / 9');
A2A2_k1_integral += 1/4 * A2A2_T9_term;


print('Result for A2A2: ')
print(A2A2_k1_integral);

# Save result in text file:
f = open('A2A2_k1_integral.txt', 'w');
f.write(str(A2A2_k1_integral));
f.close();
