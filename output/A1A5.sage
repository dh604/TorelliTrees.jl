from admcycles import *

def lambdaclass(d,g,n):
    R = TautologicalRing(g,n,moduli='ct')
    return R.lambdaclass(d)

A1A5_summands = [];
R = TautologicalRing(6,0,moduli='ct')

A1A5_graph = StableGraph([1, 5], [[1], [2]], [(1, 2)]);
A1A5_T1_terms = [];

A1A5_T1_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 1), psiclass(1, 5, 1)^4]));
print('Computed term 1 / 5 of tree 1 / 24');
A1A5_T1_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 1), psiclass(1, 5, 1)^3*lambdaclass(1, 5, 1)]));
print('Computed term 2 / 5 of tree 1 / 24');
A1A5_T1_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 1), psiclass(1, 5, 1)^2*lambdaclass(2, 5, 1)]));
print('Computed term 3 / 5 of tree 1 / 24');
A1A5_T1_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 1), psiclass(1, 5, 1)*lambdaclass(3, 5, 1)]));
print('Computed term 4 / 5 of tree 1 / 24');
A1A5_T1_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 1), lambdaclass(4, 5, 1)]));
print('Computed term 5 / 5 of tree 1 / 24');
A1A5_summands.append( 1/1 * R.sum(A1A5_T1_terms));

A1A5_graph = StableGraph([1, 4, 1], [[1, 2], [3], [4]], [(1, 3), (2, 4)]);
A1A5_T2_terms = [];

A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([3 * psiclass(1, 1, 2), psiclass(1, 4, 1)^2, fundclass(1, 1)]));
print('Computed term 1 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 2), psiclass(1, 4, 1)^2, fundclass(1, 1)]));
print('Computed term 2 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 4, 1)^3, fundclass(1, 1)]));
print('Computed term 3 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([-5 * lambdaclass(1, 1, 2), psiclass(1, 4, 1)^2, fundclass(1, 1)]));
print('Computed term 4 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([-2 * psiclass(1, 1, 2), psiclass(1, 4, 1)*lambdaclass(1, 4, 1), fundclass(1, 1)]));
print('Computed term 5 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(2, 1, 2), psiclass(1, 4, 1)*lambdaclass(1, 4, 1), fundclass(1, 1)]));
print('Computed term 6 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 2), psiclass(1, 4, 1)^2*lambdaclass(1, 4, 1), fundclass(1, 1)]));
print('Computed term 7 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([4 * lambdaclass(1, 1, 2), psiclass(1, 4, 1)*lambdaclass(1, 4, 1), fundclass(1, 1)]));
print('Computed term 8 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 2), lambdaclass(2, 4, 1), fundclass(1, 1)]));
print('Computed term 9 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 2), lambdaclass(2, 4, 1), fundclass(1, 1)]));
print('Computed term 10 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 4, 1)*lambdaclass(2, 4, 1), fundclass(1, 1)]));
print('Computed term 11 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([-3 * lambdaclass(1, 1, 2), lambdaclass(2, 4, 1), fundclass(1, 1)]));
print('Computed term 12 / 13 of tree 2 / 24');
A1A5_T2_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 2), lambdaclass(3, 4, 1), fundclass(1, 1)]));
print('Computed term 13 / 13 of tree 2 / 24');
A1A5_summands.append( 1/1 * R.sum(A1A5_T2_terms));

A1A5_graph = StableGraph([1, 3, 2], [[1, 2], [3], [4]], [(1, 3), (2, 4)]);
A1A5_T3_terms = [];

A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([3 * psiclass(1, 1, 2), psiclass(1, 3, 1)^2, fundclass(2, 1)]));
print('Computed term 1 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 2), psiclass(1, 3, 1)^2, fundclass(2, 1)]));
print('Computed term 2 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 3, 1)^3, fundclass(2, 1)]));
print('Computed term 3 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([2 * psiclass(1, 1, 2), psiclass(1, 3, 1), psiclass(1, 2, 1)]));
print('Computed term 4 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([2 * psiclass(2, 1, 2), psiclass(1, 3, 1), psiclass(1, 2, 1)]));
print('Computed term 5 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 3, 1)^2, psiclass(1, 2, 1)]));
print('Computed term 6 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 2), fundclass(3, 1), psiclass(1, 2, 1)^2]));
print('Computed term 7 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([3 * psiclass(2, 1, 2), fundclass(3, 1), psiclass(1, 2, 1)^2]));
print('Computed term 8 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 3, 1), psiclass(1, 2, 1)^2]));
print('Computed term 9 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-5 * lambdaclass(1, 1, 2), psiclass(1, 3, 1)^2, fundclass(2, 1)]));
print('Computed term 10 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-5 * lambdaclass(1, 1, 2), psiclass(1, 3, 1), psiclass(1, 2, 1)]));
print('Computed term 11 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-5 * lambdaclass(1, 1, 2), fundclass(3, 1), psiclass(1, 2, 1)^2]));
print('Computed term 12 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-2 * psiclass(1, 1, 2), psiclass(1, 3, 1)*lambdaclass(1, 3, 1), fundclass(2, 1)]));
print('Computed term 13 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(2, 1, 2), psiclass(1, 3, 1)*lambdaclass(1, 3, 1), fundclass(2, 1)]));
print('Computed term 14 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 2), psiclass(1, 3, 1)^2*lambdaclass(1, 3, 1), fundclass(2, 1)]));
print('Computed term 15 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(1, 1, 2), lambdaclass(1, 3, 1), psiclass(1, 2, 1)]));
print('Computed term 16 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-2 * psiclass(2, 1, 2), lambdaclass(1, 3, 1), psiclass(1, 2, 1)]));
print('Computed term 17 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 2), psiclass(1, 3, 1)*lambdaclass(1, 3, 1), psiclass(1, 2, 1)]));
print('Computed term 18 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 2), lambdaclass(1, 3, 1), psiclass(1, 2, 1)^2]));
print('Computed term 19 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([4 * lambdaclass(1, 1, 2), psiclass(1, 3, 1)*lambdaclass(1, 3, 1), fundclass(2, 1)]));
print('Computed term 20 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([4 * lambdaclass(1, 1, 2), lambdaclass(1, 3, 1), psiclass(1, 2, 1)]));
print('Computed term 21 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 2), lambdaclass(2, 3, 1), fundclass(2, 1)]));
print('Computed term 22 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 2), lambdaclass(2, 3, 1), fundclass(2, 1)]));
print('Computed term 23 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 3, 1)*lambdaclass(2, 3, 1), fundclass(2, 1)]));
print('Computed term 24 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 2), lambdaclass(2, 3, 1), psiclass(1, 2, 1)]));
print('Computed term 25 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-3 * lambdaclass(1, 1, 2), lambdaclass(2, 3, 1), fundclass(2, 1)]));
print('Computed term 26 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 2), lambdaclass(3, 3, 1), fundclass(2, 1)]));
print('Computed term 27 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-2 * psiclass(1, 1, 2), psiclass(1, 3, 1), lambdaclass(1, 2, 1)]));
print('Computed term 28 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(2, 1, 2), psiclass(1, 3, 1), lambdaclass(1, 2, 1)]));
print('Computed term 29 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 2), psiclass(1, 3, 1)^2, lambdaclass(1, 2, 1)]));
print('Computed term 30 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(1, 1, 2), fundclass(3, 1), psiclass(1, 2, 1)*lambdaclass(1, 2, 1)]));
print('Computed term 31 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-2 * psiclass(2, 1, 2), fundclass(3, 1), psiclass(1, 2, 1)*lambdaclass(1, 2, 1)]));
print('Computed term 32 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 2), psiclass(1, 3, 1), psiclass(1, 2, 1)*lambdaclass(1, 2, 1)]));
print('Computed term 33 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([4 * lambdaclass(1, 1, 2), psiclass(1, 3, 1), lambdaclass(1, 2, 1)]));
print('Computed term 34 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([4 * lambdaclass(1, 1, 2), fundclass(3, 1), psiclass(1, 2, 1)*lambdaclass(1, 2, 1)]));
print('Computed term 35 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 2), lambdaclass(1, 3, 1), lambdaclass(1, 2, 1)]));
print('Computed term 36 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 2), lambdaclass(1, 3, 1), lambdaclass(1, 2, 1)]));
print('Computed term 37 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 3, 1)*lambdaclass(1, 3, 1), lambdaclass(1, 2, 1)]));
print('Computed term 38 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 2), lambdaclass(1, 3, 1), psiclass(1, 2, 1)*lambdaclass(1, 2, 1)]));
print('Computed term 39 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-3 * lambdaclass(1, 1, 2), lambdaclass(1, 3, 1), lambdaclass(1, 2, 1)]));
print('Computed term 40 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 2), lambdaclass(2, 3, 1), lambdaclass(1, 2, 1)]));
print('Computed term 41 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 2), fundclass(3, 1), lambdaclass(2, 2, 1)]));
print('Computed term 42 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 2), fundclass(3, 1), lambdaclass(2, 2, 1)]));
print('Computed term 43 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 3, 1), lambdaclass(2, 2, 1)]));
print('Computed term 44 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-3 * lambdaclass(1, 1, 2), fundclass(3, 1), lambdaclass(2, 2, 1)]));
print('Computed term 45 / 46 of tree 3 / 24');
A1A5_T3_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 2), lambdaclass(1, 3, 1), lambdaclass(2, 2, 1)]));
print('Computed term 46 / 46 of tree 3 / 24');
A1A5_summands.append( 1/1 * R.sum(A1A5_T3_terms));

A1A5_graph = StableGraph([0, 1, 4, 1], [[1, 2, 3], [4], [5], [6]], [(1, 4), (2, 5), (3, 6)]);
A1A5_T4_terms = [];

A1A5_T4_terms.append(A1A5_graph.boundary_pushforward([-5 * fundclass(0, 3), fundclass(1, 1), psiclass(1, 4, 1)^2, fundclass(1, 1)]));
print('Computed term 1 / 3 of tree 4 / 24');
A1A5_T4_terms.append(A1A5_graph.boundary_pushforward([4 * fundclass(0, 3), fundclass(1, 1), psiclass(1, 4, 1)*lambdaclass(1, 4, 1), fundclass(1, 1)]));
print('Computed term 2 / 3 of tree 4 / 24');
A1A5_T4_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), fundclass(1, 1), lambdaclass(2, 4, 1), fundclass(1, 1)]));
print('Computed term 3 / 3 of tree 4 / 24');
A1A5_summands.append( 1/1 * R.sum(A1A5_T4_terms));

A1A5_graph = StableGraph([0, 1, 3, 2], [[1, 2, 3], [4], [5], [6]], [(1, 4), (2, 5), (3, 6)]);
A1A5_T5_terms = [];

A1A5_T5_terms.append(A1A5_graph.boundary_pushforward([-5 * fundclass(0, 3), fundclass(1, 1), psiclass(1, 3, 1)^2, fundclass(2, 1)]));
print('Computed term 1 / 10 of tree 5 / 24');
A1A5_T5_terms.append(A1A5_graph.boundary_pushforward([-5 * fundclass(0, 3), fundclass(1, 1), psiclass(1, 3, 1), psiclass(1, 2, 1)]));
print('Computed term 2 / 10 of tree 5 / 24');
A1A5_T5_terms.append(A1A5_graph.boundary_pushforward([-5 * fundclass(0, 3), fundclass(1, 1), fundclass(3, 1), psiclass(1, 2, 1)^2]));
print('Computed term 3 / 10 of tree 5 / 24');
A1A5_T5_terms.append(A1A5_graph.boundary_pushforward([4 * fundclass(0, 3), fundclass(1, 1), psiclass(1, 3, 1)*lambdaclass(1, 3, 1), fundclass(2, 1)]));
print('Computed term 4 / 10 of tree 5 / 24');
A1A5_T5_terms.append(A1A5_graph.boundary_pushforward([4 * fundclass(0, 3), fundclass(1, 1), lambdaclass(1, 3, 1), psiclass(1, 2, 1)]));
print('Computed term 5 / 10 of tree 5 / 24');
A1A5_T5_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), fundclass(1, 1), lambdaclass(2, 3, 1), fundclass(2, 1)]));
print('Computed term 6 / 10 of tree 5 / 24');
A1A5_T5_terms.append(A1A5_graph.boundary_pushforward([4 * fundclass(0, 3), fundclass(1, 1), psiclass(1, 3, 1), lambdaclass(1, 2, 1)]));
print('Computed term 7 / 10 of tree 5 / 24');
A1A5_T5_terms.append(A1A5_graph.boundary_pushforward([4 * fundclass(0, 3), fundclass(1, 1), fundclass(3, 1), psiclass(1, 2, 1)*lambdaclass(1, 2, 1)]));
print('Computed term 8 / 10 of tree 5 / 24');
A1A5_T5_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), fundclass(1, 1), lambdaclass(1, 3, 1), lambdaclass(1, 2, 1)]));
print('Computed term 9 / 10 of tree 5 / 24');
A1A5_T5_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), fundclass(1, 1), fundclass(3, 1), lambdaclass(2, 2, 1)]));
print('Computed term 10 / 10 of tree 5 / 24');
A1A5_summands.append( 1/1 * R.sum(A1A5_T5_terms));

A1A5_graph = StableGraph([1, 3, 1, 1], [[1, 2, 3], [4], [5], [6]], [(1, 4), (2, 5), (3, 6)]);
A1A5_T6_terms = [];

A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 3)^2, fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 3)*psiclass(2, 1, 3), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 2 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 3)^2, fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 3 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 3)*psiclass(3, 1, 3), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 4 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 3)*psiclass(3, 1, 3), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 5 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(3, 1, 3)^2, fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 6 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([2 * psiclass(1, 1, 3), psiclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 7 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 3), psiclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 8 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(3, 1, 3), psiclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 9 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 3), psiclass(1, 3, 1)^2, fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 10 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(1, 1, 3)*lambdaclass(1, 1, 3), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 11 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(2, 1, 3)*lambdaclass(1, 1, 3), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 12 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(3, 1, 3)*lambdaclass(1, 1, 3), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 13 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([-5 * lambdaclass(1, 1, 3), psiclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 14 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([10 * lambdaclass(1, 1, 3)^2, fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 15 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(1, 1, 3), lambdaclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 16 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(2, 1, 3), lambdaclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 17 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(3, 1, 3), lambdaclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 18 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 3), psiclass(1, 3, 1)*lambdaclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 19 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([4 * lambdaclass(1, 1, 3), lambdaclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 20 / 21 of tree 6 / 24');
A1A5_T6_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 3), lambdaclass(2, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 21 / 21 of tree 6 / 24');
A1A5_summands.append( 1/2 * R.sum(A1A5_T6_terms));

A1A5_graph = StableGraph([1, 2, 2, 1], [[1, 2, 3], [4], [5], [6]], [(1, 4), (2, 5), (3, 6)]);
A1A5_T7_terms = [];

A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 3)^2, fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 1 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 3)*psiclass(2, 1, 3), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 2 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 3)^2, fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 3 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 3)*psiclass(3, 1, 3), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 4 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 3)*psiclass(3, 1, 3), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 5 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(3, 1, 3)^2, fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 6 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([2 * psiclass(1, 1, 3), psiclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 7 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 3), psiclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 8 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(3, 1, 3), psiclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 9 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 3), psiclass(1, 2, 1)^2, fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 10 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 3), fundclass(2, 1), psiclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 11 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([2 * psiclass(2, 1, 3), fundclass(2, 1), psiclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 12 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(3, 1, 3), fundclass(2, 1), psiclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 13 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 3), psiclass(1, 2, 1), psiclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 14 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 3), fundclass(2, 1), psiclass(1, 2, 1)^2, fundclass(1, 1)]));
print('Computed term 15 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(1, 1, 3)*lambdaclass(1, 1, 3), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 16 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(2, 1, 3)*lambdaclass(1, 1, 3), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 17 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(3, 1, 3)*lambdaclass(1, 1, 3), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 18 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-5 * lambdaclass(1, 1, 3), psiclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 19 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-5 * lambdaclass(1, 1, 3), fundclass(2, 1), psiclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 20 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([10 * lambdaclass(1, 1, 3)^2, fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 21 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(1, 1, 3), lambdaclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 22 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(2, 1, 3), lambdaclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 23 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(3, 1, 3), lambdaclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 24 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 3), psiclass(1, 2, 1)*lambdaclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 25 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 3), lambdaclass(1, 2, 1), psiclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 26 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([4 * lambdaclass(1, 1, 3), lambdaclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 27 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 3), lambdaclass(2, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 28 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(1, 1, 3), fundclass(2, 1), lambdaclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 29 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(2, 1, 3), fundclass(2, 1), lambdaclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 30 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-1 * psiclass(3, 1, 3), fundclass(2, 1), lambdaclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 31 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 3), psiclass(1, 2, 1), lambdaclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 32 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 3), fundclass(2, 1), psiclass(1, 2, 1)*lambdaclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 33 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([4 * lambdaclass(1, 1, 3), fundclass(2, 1), lambdaclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 34 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 3), lambdaclass(1, 2, 1), lambdaclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 35 / 36 of tree 7 / 24');
A1A5_T7_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 3), fundclass(2, 1), lambdaclass(2, 2, 1), fundclass(1, 1)]));
print('Computed term 36 / 36 of tree 7 / 24');
A1A5_summands.append( 1/2 * R.sum(A1A5_T7_terms));

A1A5_graph = StableGraph([0, 1, 3, 1, 1], [[1, 2, 3], [4, 5], [6], [7], [8]], [(1, 4), (2, 7), (3, 8), (5, 6)]);
A1A5_T8_terms = [];

A1A5_T8_terms.append(A1A5_graph.boundary_pushforward([-6 * fundclass(0, 3), psiclass(1, 1, 2), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 5 of tree 8 / 24');
A1A5_T8_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), psiclass(2, 1, 2), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 2 / 5 of tree 8 / 24');
A1A5_T8_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), fundclass(1, 2), psiclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 3 / 5 of tree 8 / 24');
A1A5_T8_terms.append(A1A5_graph.boundary_pushforward([15 * fundclass(0, 3), lambdaclass(1, 1, 2), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 4 / 5 of tree 8 / 24');
A1A5_T8_terms.append(A1A5_graph.boundary_pushforward([3 * fundclass(0, 3), fundclass(1, 2), lambdaclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 5 / 5 of tree 8 / 24');
A1A5_summands.append( 1/2 * R.sum(A1A5_T8_terms));

A1A5_graph = StableGraph([0, 1, 2, 2, 1], [[1, 2, 3], [4, 5], [6], [7], [8]], [(1, 4), (2, 7), (3, 8), (5, 6)]);
A1A5_T9_terms = [];

A1A5_T9_terms.append(A1A5_graph.boundary_pushforward([-6 * fundclass(0, 3), psiclass(1, 1, 2), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 1 / 7 of tree 9 / 24');
A1A5_T9_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), psiclass(2, 1, 2), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 2 / 7 of tree 9 / 24');
A1A5_T9_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), fundclass(1, 2), psiclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 3 / 7 of tree 9 / 24');
A1A5_T9_terms.append(A1A5_graph.boundary_pushforward([-4 * fundclass(0, 3), fundclass(1, 2), fundclass(2, 1), psiclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 4 / 7 of tree 9 / 24');
A1A5_T9_terms.append(A1A5_graph.boundary_pushforward([15 * fundclass(0, 3), lambdaclass(1, 1, 2), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 5 / 7 of tree 9 / 24');
A1A5_T9_terms.append(A1A5_graph.boundary_pushforward([3 * fundclass(0, 3), fundclass(1, 2), lambdaclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 6 / 7 of tree 9 / 24');
A1A5_T9_terms.append(A1A5_graph.boundary_pushforward([3 * fundclass(0, 3), fundclass(1, 2), fundclass(2, 1), lambdaclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 7 / 7 of tree 9 / 24');
A1A5_summands.append( 1/1 * R.sum(A1A5_T9_terms));

A1A5_graph = StableGraph([0, 1, 1, 3, 1], [[1, 2, 3], [4, 5], [6], [7], [8]], [(1, 4), (2, 7), (3, 8), (5, 6)]);
A1A5_T10_terms = [];

A1A5_T10_terms.append(A1A5_graph.boundary_pushforward([-6 * fundclass(0, 3), psiclass(1, 1, 2), fundclass(1, 1), fundclass(3, 1), fundclass(1, 1)]));
print('Computed term 1 / 5 of tree 10 / 24');
A1A5_T10_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), psiclass(2, 1, 2), fundclass(1, 1), fundclass(3, 1), fundclass(1, 1)]));
print('Computed term 2 / 5 of tree 10 / 24');
A1A5_T10_terms.append(A1A5_graph.boundary_pushforward([-4 * fundclass(0, 3), fundclass(1, 2), fundclass(1, 1), psiclass(1, 3, 1), fundclass(1, 1)]));
print('Computed term 3 / 5 of tree 10 / 24');
A1A5_T10_terms.append(A1A5_graph.boundary_pushforward([15 * fundclass(0, 3), lambdaclass(1, 1, 2), fundclass(1, 1), fundclass(3, 1), fundclass(1, 1)]));
print('Computed term 4 / 5 of tree 10 / 24');
A1A5_T10_terms.append(A1A5_graph.boundary_pushforward([3 * fundclass(0, 3), fundclass(1, 2), fundclass(1, 1), lambdaclass(1, 3, 1), fundclass(1, 1)]));
print('Computed term 5 / 5 of tree 10 / 24');
A1A5_summands.append( 1/1 * R.sum(A1A5_T10_terms));

A1A5_graph = StableGraph([0, 1, 1, 2, 2], [[1, 2, 3], [4, 5], [6], [7], [8]], [(1, 4), (2, 7), (3, 8), (5, 6)]);
A1A5_T11_terms = [];

A1A5_T11_terms.append(A1A5_graph.boundary_pushforward([-6 * fundclass(0, 3), psiclass(1, 1, 2), fundclass(1, 1), fundclass(2, 1), fundclass(2, 1)]));
print('Computed term 1 / 7 of tree 11 / 24');
A1A5_T11_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), psiclass(2, 1, 2), fundclass(1, 1), fundclass(2, 1), fundclass(2, 1)]));
print('Computed term 2 / 7 of tree 11 / 24');
A1A5_T11_terms.append(A1A5_graph.boundary_pushforward([-4 * fundclass(0, 3), fundclass(1, 2), fundclass(1, 1), psiclass(1, 2, 1), fundclass(2, 1)]));
print('Computed term 3 / 7 of tree 11 / 24');
A1A5_T11_terms.append(A1A5_graph.boundary_pushforward([-4 * fundclass(0, 3), fundclass(1, 2), fundclass(1, 1), fundclass(2, 1), psiclass(1, 2, 1)]));
print('Computed term 4 / 7 of tree 11 / 24');
A1A5_T11_terms.append(A1A5_graph.boundary_pushforward([15 * fundclass(0, 3), lambdaclass(1, 1, 2), fundclass(1, 1), fundclass(2, 1), fundclass(2, 1)]));
print('Computed term 5 / 7 of tree 11 / 24');
A1A5_T11_terms.append(A1A5_graph.boundary_pushforward([3 * fundclass(0, 3), fundclass(1, 2), fundclass(1, 1), lambdaclass(1, 2, 1), fundclass(2, 1)]));
print('Computed term 6 / 7 of tree 11 / 24');
A1A5_T11_terms.append(A1A5_graph.boundary_pushforward([3 * fundclass(0, 3), fundclass(1, 2), fundclass(1, 1), fundclass(2, 1), lambdaclass(1, 2, 1)]));
print('Computed term 7 / 7 of tree 11 / 24');
A1A5_summands.append( 1/2 * R.sum(A1A5_T11_terms));

A1A5_graph = StableGraph([0, 1, 3, 1, 1], [[1, 2, 3, 4], [5], [6], [7], [8]], [(1, 5), (2, 6), (3, 7), (4, 8)]);
A1A5_T12_terms = [];

A1A5_T12_terms.append(A1A5_graph.boundary_pushforward([-10 * psiclass(1, 0, 4), fundclass(1, 1), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 6 of tree 12 / 24');
A1A5_T12_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(2, 0, 4), fundclass(1, 1), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 2 / 6 of tree 12 / 24');
A1A5_T12_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(3, 0, 4), fundclass(1, 1), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 3 / 6 of tree 12 / 24');
A1A5_T12_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(4, 0, 4), fundclass(1, 1), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 4 / 6 of tree 12 / 24');
A1A5_T12_terms.append(A1A5_graph.boundary_pushforward([-5 * fundclass(0, 4), fundclass(1, 1), psiclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 5 / 6 of tree 12 / 24');
A1A5_T12_terms.append(A1A5_graph.boundary_pushforward([4 * fundclass(0, 4), fundclass(1, 1), lambdaclass(1, 3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 6 / 6 of tree 12 / 24');
A1A5_summands.append( 1/2 * R.sum(A1A5_T12_terms));

A1A5_graph = StableGraph([0, 1, 2, 2, 1], [[1, 2, 3, 4], [5], [6], [7], [8]], [(1, 5), (2, 6), (3, 7), (4, 8)]);
A1A5_T13_terms = [];

A1A5_T13_terms.append(A1A5_graph.boundary_pushforward([-10 * psiclass(1, 0, 4), fundclass(1, 1), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 1 / 8 of tree 13 / 24');
A1A5_T13_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(2, 0, 4), fundclass(1, 1), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 2 / 8 of tree 13 / 24');
A1A5_T13_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(3, 0, 4), fundclass(1, 1), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 3 / 8 of tree 13 / 24');
A1A5_T13_terms.append(A1A5_graph.boundary_pushforward([-5 * psiclass(4, 0, 4), fundclass(1, 1), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 4 / 8 of tree 13 / 24');
A1A5_T13_terms.append(A1A5_graph.boundary_pushforward([-5 * fundclass(0, 4), fundclass(1, 1), psiclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 5 / 8 of tree 13 / 24');
A1A5_T13_terms.append(A1A5_graph.boundary_pushforward([-5 * fundclass(0, 4), fundclass(1, 1), fundclass(2, 1), psiclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 6 / 8 of tree 13 / 24');
A1A5_T13_terms.append(A1A5_graph.boundary_pushforward([4 * fundclass(0, 4), fundclass(1, 1), lambdaclass(1, 2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 7 / 8 of tree 13 / 24');
A1A5_T13_terms.append(A1A5_graph.boundary_pushforward([4 * fundclass(0, 4), fundclass(1, 1), fundclass(2, 1), lambdaclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 8 / 8 of tree 13 / 24');
A1A5_summands.append( 1/2 * R.sum(A1A5_T13_terms));

A1A5_graph = StableGraph([1, 2, 1, 1, 1], [[1, 2, 3, 4], [5], [6], [7], [8]], [(1, 5), (2, 6), (3, 7), (4, 8)]);
A1A5_T14_terms = [];

A1A5_T14_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(1, 1, 4), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 7 of tree 14 / 24');
A1A5_T14_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(2, 1, 4), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 2 / 7 of tree 14 / 24');
A1A5_T14_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(3, 1, 4), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 3 / 7 of tree 14 / 24');
A1A5_T14_terms.append(A1A5_graph.boundary_pushforward([1 * psiclass(4, 1, 4), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 4 / 7 of tree 14 / 24');
A1A5_T14_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 4), psiclass(1, 2, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 5 / 7 of tree 14 / 24');
A1A5_T14_terms.append(A1A5_graph.boundary_pushforward([-5 * lambdaclass(1, 1, 4), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 6 / 7 of tree 14 / 24');
A1A5_T14_terms.append(A1A5_graph.boundary_pushforward([-1 * fundclass(1, 4), lambdaclass(1, 2, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 7 / 7 of tree 14 / 24');
A1A5_summands.append( 1/6 * R.sum(A1A5_T14_terms));

A1A5_graph = StableGraph([0, 0, 1, 3, 1, 1], [[1, 2, 3], [4, 5, 6], [7], [8], [9], [10]], [(1, 4), (2, 9), (3, 10), (5, 7), (6, 8)]);
A1A5_T15_terms = [];

A1A5_T15_terms.append(A1A5_graph.boundary_pushforward([15 * fundclass(0, 3), fundclass(0, 3), fundclass(1, 1), fundclass(3, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 15 / 24');
A1A5_summands.append( 1/2 * R.sum(A1A5_T15_terms));

A1A5_graph = StableGraph([0, 0, 1, 2, 2, 1], [[1, 2, 3], [4, 5, 6], [7], [8], [9], [10]], [(1, 4), (2, 9), (3, 10), (5, 7), (6, 8)]);
A1A5_T16_terms = [];

A1A5_T16_terms.append(A1A5_graph.boundary_pushforward([15 * fundclass(0, 3), fundclass(0, 3), fundclass(1, 1), fundclass(2, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 16 / 24');
A1A5_summands.append( 1/1 * R.sum(A1A5_T16_terms));

A1A5_graph = StableGraph([0, 0, 1, 1, 3, 1], [[1, 2, 3], [4, 5, 6], [7], [8], [9], [10]], [(1, 4), (2, 9), (3, 10), (5, 7), (6, 8)]);
A1A5_T17_terms = [];

A1A5_T17_terms.append(A1A5_graph.boundary_pushforward([15 * fundclass(0, 3), fundclass(0, 3), fundclass(1, 1), fundclass(1, 1), fundclass(3, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 17 / 24');
A1A5_summands.append( 1/1 * R.sum(A1A5_T17_terms));

A1A5_graph = StableGraph([0, 0, 1, 1, 2, 2], [[1, 2, 3], [4, 5, 6], [7], [8], [9], [10]], [(1, 4), (2, 9), (3, 10), (5, 7), (6, 8)]);
A1A5_T18_terms = [];

A1A5_T18_terms.append(A1A5_graph.boundary_pushforward([15 * fundclass(0, 3), fundclass(0, 3), fundclass(1, 1), fundclass(1, 1), fundclass(2, 1), fundclass(2, 1)]));
print('Computed term 1 / 1 of tree 18 / 24');
A1A5_summands.append( 1/2 * R.sum(A1A5_T18_terms));

A1A5_graph = StableGraph([0, 1, 2, 1, 1, 1], [[1, 2, 3], [4, 5, 6], [7], [8], [9], [10]], [(1, 4), (2, 9), (3, 10), (5, 7), (6, 8)]);
A1A5_T19_terms = [];

A1A5_T19_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), fundclass(1, 3), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 19 / 24');
A1A5_summands.append( 1/2 * R.sum(A1A5_T19_terms));

A1A5_graph = StableGraph([0, 1, 1, 1, 2, 1], [[1, 2, 3], [4, 5, 6], [7], [8], [9], [10]], [(1, 4), (2, 9), (3, 10), (5, 7), (6, 8)]);
A1A5_T20_terms = [];

A1A5_T20_terms.append(A1A5_graph.boundary_pushforward([-3 * fundclass(0, 3), fundclass(1, 3), fundclass(1, 1), fundclass(1, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 20 / 24');
A1A5_summands.append( 1/2 * R.sum(A1A5_T20_terms));

A1A5_graph = StableGraph([0, 1, 2, 1, 1, 1], [[1, 2, 3, 4], [5, 6], [7], [8], [9], [10]], [(1, 5), (2, 8), (3, 9), (4, 10), (6, 7)]);
A1A5_T21_terms = [];

A1A5_T21_terms.append(A1A5_graph.boundary_pushforward([-4 * fundclass(0, 4), fundclass(1, 2), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 21 / 24');
A1A5_summands.append( 1/6 * R.sum(A1A5_T21_terms));

A1A5_graph = StableGraph([0, 1, 1, 2, 1, 1], [[1, 2, 3, 4], [5, 6], [7], [8], [9], [10]], [(1, 5), (2, 8), (3, 9), (4, 10), (6, 7)]);
A1A5_T22_terms = [];

A1A5_T22_terms.append(A1A5_graph.boundary_pushforward([-4 * fundclass(0, 4), fundclass(1, 2), fundclass(1, 1), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 22 / 24');
A1A5_summands.append( 1/2 * R.sum(A1A5_T22_terms));

A1A5_graph = StableGraph([0, 1, 2, 1, 1, 1], [[1, 2, 3, 4, 5], [6], [7], [8], [9], [10]], [(1, 6), (2, 7), (3, 8), (4, 9), (5, 10)]);
A1A5_T23_terms = [];

A1A5_T23_terms.append(A1A5_graph.boundary_pushforward([-5 * fundclass(0, 5), fundclass(1, 1), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 23 / 24');
A1A5_summands.append( 1/6 * R.sum(A1A5_T23_terms));

A1A5_graph = StableGraph([1, 1, 1, 1, 1, 1], [[1, 2, 3, 4, 5], [6], [7], [8], [9], [10]], [(1, 6), (2, 7), (3, 8), (4, 9), (5, 10)]);
A1A5_T24_terms = [];

A1A5_T24_terms.append(A1A5_graph.boundary_pushforward([1 * fundclass(1, 5), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 24 / 24');
A1A5_summands.append( 1/120 * R.sum(A1A5_T24_terms));

Torelli_pullback = R.sum(A1A5_summands);
