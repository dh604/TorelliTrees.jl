from admcycles import *

def lambdaclass(d,g,n):
    R = TautologicalRing(g,n,moduli='ct')
    return R.lambdaclass(d)

A1A4_summands = [];
R = TautologicalRing(5,0,moduli='ct')

A1A4_graph = StableGraph([1, 4], [[1], [2]], [(1, 2)]);
A1A4_T1_terms = [];

A1A4_T1_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 1), psiclass(1, 4, 1)^3]));
print('Computed term 1 / 4 of tree 1 / 10');
A1A4_T1_terms.append(A1A4_graph.boundary_pushforward([-1 * fundclass(1, 1), psiclass(1, 4, 1)^2*lambdaclass(1, 4, 1)]));
print('Computed term 2 / 4 of tree 1 / 10');
A1A4_T1_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 1), psiclass(1, 4, 1)*lambdaclass(2, 4, 1)]));
print('Computed term 3 / 4 of tree 1 / 10');
A1A4_T1_terms.append(A1A4_graph.boundary_pushforward([-1 * fundclass(1, 1), lambdaclass(3, 4, 1)]));
print('Computed term 4 / 4 of tree 1 / 10');
A1A4_summands.append( 1/1 * R.sum(A1A4_T1_terms));

A1A4_graph = StableGraph([1, 3, 1], [[1, 2], [3], [4]], [(1, 3), (2, 4)]);
A1A4_T2_terms = [];

A1A4_T2_terms.append(A1A4_graph.boundary_pushforward([2 * psiclass(1, 1, 2), psiclass(1, 3, 1), fundclass(1, 1)]));
print('Computed term 1 / 9 of tree 2 / 10');
A1A4_T2_terms.append(A1A4_graph.boundary_pushforward([1 * psiclass(2, 1, 2), psiclass(1, 3, 1), fundclass(1, 1)]));
print('Computed term 2 / 9 of tree 2 / 10');
A1A4_T2_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 3, 1)^2, fundclass(1, 1)]));
print('Computed term 3 / 9 of tree 2 / 10');
A1A4_T2_terms.append(A1A4_graph.boundary_pushforward([-4 * lambdaclass(1, 1, 2), psiclass(1, 3, 1), fundclass(1, 1)]));
print('Computed term 4 / 9 of tree 2 / 10');
A1A4_T2_terms.append(A1A4_graph.boundary_pushforward([-1 * psiclass(1, 1, 2), lambdaclass(1, 3, 1), fundclass(1, 1)]));
print('Computed term 5 / 9 of tree 2 / 10');
A1A4_T2_terms.append(A1A4_graph.boundary_pushforward([-1 * psiclass(2, 1, 2), lambdaclass(1, 3, 1), fundclass(1, 1)]));
print('Computed term 6 / 9 of tree 2 / 10');
A1A4_T2_terms.append(A1A4_graph.boundary_pushforward([-1 * fundclass(1, 2), psiclass(1, 3, 1)*lambdaclass(1, 3, 1), fundclass(1, 1)]));
print('Computed term 7 / 9 of tree 2 / 10');
A1A4_T2_terms.append(A1A4_graph.boundary_pushforward([3 * lambdaclass(1, 1, 2), lambdaclass(1, 3, 1), fundclass(1, 1)]));
print('Computed term 8 / 9 of tree 2 / 10');
A1A4_T2_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 2), lambdaclass(2, 3, 1), fundclass(1, 1)]));
print('Computed term 9 / 9 of tree 2 / 10');
A1A4_summands.append( 1/1 * R.sum(A1A4_T2_terms));

A1A4_graph = StableGraph([1, 2, 2], [[1, 2], [3], [4]], [(1, 3), (2, 4)]);
A1A4_T3_terms = [];

A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([2 * psiclass(1, 1, 2), psiclass(1, 2, 1), fundclass(2, 1)]));
print('Computed term 1 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([1 * psiclass(2, 1, 2), psiclass(1, 2, 1), fundclass(2, 1)]));
print('Computed term 2 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 2, 1)^2, fundclass(2, 1)]));
print('Computed term 3 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([1 * psiclass(1, 1, 2), fundclass(2, 1), psiclass(1, 2, 1)]));
print('Computed term 4 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([2 * psiclass(2, 1, 2), fundclass(2, 1), psiclass(1, 2, 1)]));
print('Computed term 5 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 2, 1), psiclass(1, 2, 1)]));
print('Computed term 6 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 2), fundclass(2, 1), psiclass(1, 2, 1)^2]));
print('Computed term 7 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([-4 * lambdaclass(1, 1, 2), psiclass(1, 2, 1), fundclass(2, 1)]));
print('Computed term 8 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([-4 * lambdaclass(1, 1, 2), fundclass(2, 1), psiclass(1, 2, 1)]));
print('Computed term 9 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([-1 * psiclass(1, 1, 2), lambdaclass(1, 2, 1), fundclass(2, 1)]));
print('Computed term 10 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([-1 * psiclass(2, 1, 2), lambdaclass(1, 2, 1), fundclass(2, 1)]));
print('Computed term 11 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([-1 * fundclass(1, 2), psiclass(1, 2, 1)*lambdaclass(1, 2, 1), fundclass(2, 1)]));
print('Computed term 12 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([-1 * fundclass(1, 2), lambdaclass(1, 2, 1), psiclass(1, 2, 1)]));
print('Computed term 13 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([3 * lambdaclass(1, 1, 2), lambdaclass(1, 2, 1), fundclass(2, 1)]));
print('Computed term 14 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 2), lambdaclass(2, 2, 1), fundclass(2, 1)]));
print('Computed term 15 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([-1 * psiclass(1, 1, 2), fundclass(2, 1), lambdaclass(1, 2, 1)]));
print('Computed term 16 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([-1 * psiclass(2, 1, 2), fundclass(2, 1), lambdaclass(1, 2, 1)]));
print('Computed term 17 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([-1 * fundclass(1, 2), psiclass(1, 2, 1), lambdaclass(1, 2, 1)]));
print('Computed term 18 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([-1 * fundclass(1, 2), fundclass(2, 1), psiclass(1, 2, 1)*lambdaclass(1, 2, 1)]));
print('Computed term 19 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([3 * lambdaclass(1, 1, 2), fundclass(2, 1), lambdaclass(1, 2, 1)]));
print('Computed term 20 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 2), lambdaclass(1, 2, 1), lambdaclass(1, 2, 1)]));
print('Computed term 21 / 22 of tree 3 / 10');
A1A4_T3_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 2), fundclass(2, 1), lambdaclass(2, 2, 1)]));
print('Computed term 22 / 22 of tree 3 / 10');
A1A4_summands.append( 1/2 * R.sum(A1A4_T3_terms));

A1A4_graph = StableGraph([0, 1, 3, 1], [[1, 2, 3], [4], [5], [6]], [(1, 4), (2, 5), (3, 6)]);
A1A4_T4_terms = [];

A1A4_T4_terms.append(A1A4_graph.boundary_pushforward([-4 * fundclass(0, 3), fundclass(1, 1), psiclass(1, 3, 1), fundclass(1, 1)]));
print('Computed term 1 / 2 of tree 4 / 10');
A1A4_T4_terms.append(A1A4_graph.boundary_pushforward([3 * fundclass(0, 3), fundclass(1, 1), lambdaclass(1, 3, 1), fundclass(1, 1)]));
print('Computed term 2 / 2 of tree 4 / 10');
A1A4_summands.append( 1/1 * R.sum(A1A4_T4_terms));

A1A4_graph = StableGraph([0, 1, 2, 2], [[1, 2, 3], [4], [5], [6]], [(1, 4), (2, 5), (3, 6)]);
A1A4_T5_terms = [];

A1A4_T5_terms.append(A1A4_graph.boundary_pushforward([-4 * fundclass(0, 3), fundclass(1, 1), psiclass(1, 2, 1), fundclass(2, 1)]));
print('Computed term 1 / 4 of tree 5 / 10');
A1A4_T5_terms.append(A1A4_graph.boundary_pushforward([-4 * fundclass(0, 3), fundclass(1, 1), fundclass(2, 1), psiclass(1, 2, 1)]));
print('Computed term 2 / 4 of tree 5 / 10');
A1A4_T5_terms.append(A1A4_graph.boundary_pushforward([3 * fundclass(0, 3), fundclass(1, 1), lambdaclass(1, 2, 1), fundclass(2, 1)]));
print('Computed term 3 / 4 of tree 5 / 10');
A1A4_T5_terms.append(A1A4_graph.boundary_pushforward([3 * fundclass(0, 3), fundclass(1, 1), fundclass(2, 1), lambdaclass(1, 2, 1)]));
print('Computed term 4 / 4 of tree 5 / 10');
A1A4_summands.append( 1/2 * R.sum(A1A4_T5_terms));

A1A4_graph = StableGraph([1, 2, 1, 1], [[1, 2, 3], [4], [5], [6]], [(1, 4), (2, 5), (3, 6)]);
A1A4_T6_terms = [];

A1A4_T6_terms.append(A1A4_graph.boundary_pushforward([1 * psiclass(1, 1, 3), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 6 of tree 6 / 10');
A1A4_T6_terms.append(A1A4_graph.boundary_pushforward([1 * psiclass(2, 1, 3), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 2 / 6 of tree 6 / 10');
A1A4_T6_terms.append(A1A4_graph.boundary_pushforward([1 * psiclass(3, 1, 3), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 3 / 6 of tree 6 / 10');
A1A4_T6_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 3), psiclass(1, 2, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 4 / 6 of tree 6 / 10');
A1A4_T6_terms.append(A1A4_graph.boundary_pushforward([-4 * lambdaclass(1, 1, 3), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 5 / 6 of tree 6 / 10');
A1A4_T6_terms.append(A1A4_graph.boundary_pushforward([-1 * fundclass(1, 3), lambdaclass(1, 2, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 6 / 6 of tree 6 / 10');
A1A4_summands.append( 1/2 * R.sum(A1A4_T6_terms));

A1A4_graph = StableGraph([0, 1, 2, 1, 1], [[1, 2, 3], [4, 5], [6], [7], [8]], [(1, 4), (2, 7), (3, 8), (5, 6)]);
A1A4_T7_terms = [];

A1A4_T7_terms.append(A1A4_graph.boundary_pushforward([-3 * fundclass(0, 3), fundclass(1, 2), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 7 / 10');
A1A4_summands.append( 1/2 * R.sum(A1A4_T7_terms));

A1A4_graph = StableGraph([0, 1, 1, 2, 1], [[1, 2, 3], [4, 5], [6], [7], [8]], [(1, 4), (2, 7), (3, 8), (5, 6)]);
A1A4_T8_terms = [];

A1A4_T8_terms.append(A1A4_graph.boundary_pushforward([-3 * fundclass(0, 3), fundclass(1, 2), fundclass(1, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 8 / 10');
A1A4_summands.append( 1/1 * R.sum(A1A4_T8_terms));

A1A4_graph = StableGraph([0, 1, 2, 1, 1], [[1, 2, 3, 4], [5], [6], [7], [8]], [(1, 5), (2, 6), (3, 7), (4, 8)]);
A1A4_T9_terms = [];

A1A4_T9_terms.append(A1A4_graph.boundary_pushforward([-4 * fundclass(0, 4), fundclass(1, 1), fundclass(2, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 9 / 10');
A1A4_summands.append( 1/2 * R.sum(A1A4_T9_terms));

A1A4_graph = StableGraph([1, 1, 1, 1, 1], [[1, 2, 3, 4], [5], [6], [7], [8]], [(1, 5), (2, 6), (3, 7), (4, 8)]);
A1A4_T10_terms = [];

A1A4_T10_terms.append(A1A4_graph.boundary_pushforward([1 * fundclass(1, 4), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 10 / 10');
A1A4_summands.append( 1/24 * R.sum(A1A4_T10_terms));

Torelli_pullback = R.sum(A1A4_summands);
