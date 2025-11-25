from admcycles import *

def lambdaclass(d,g,n):
    R = TautologicalRing(g,n,moduli='ct')
    return R.lambdaclass(d)

A1A3_summands = [];
R = TautologicalRing(4,0,moduli='ct')

A1A3_graph = StableGraph([1, 3], [[1], [2]], [(1, 2)]);
A1A3_T1_terms = [];

A1A3_T1_terms.append(A1A3_graph.boundary_pushforward([1 * fundclass(1, 1), psiclass(1, 3, 1)^2]));
print('Computed term 1 / 3 of tree 1 / 4');
A1A3_T1_terms.append(A1A3_graph.boundary_pushforward([-1 * fundclass(1, 1), psiclass(1, 3, 1)*lambdaclass(1, 3, 1)]));
print('Computed term 2 / 3 of tree 1 / 4');
A1A3_T1_terms.append(A1A3_graph.boundary_pushforward([1 * fundclass(1, 1), lambdaclass(2, 3, 1)]));
print('Computed term 3 / 3 of tree 1 / 4');
A1A3_summands.append( 1/1 * R.sum(A1A3_T1_terms));

A1A3_graph = StableGraph([1, 2, 1], [[1, 2], [3], [4]], [(1, 3), (2, 4)]);
A1A3_T2_terms = [];

A1A3_T2_terms.append(A1A3_graph.boundary_pushforward([1 * psiclass(1, 1, 2), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 1 / 5 of tree 2 / 4');
A1A3_T2_terms.append(A1A3_graph.boundary_pushforward([1 * psiclass(2, 1, 2), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 2 / 5 of tree 2 / 4');
A1A3_T2_terms.append(A1A3_graph.boundary_pushforward([1 * fundclass(1, 2), psiclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 3 / 5 of tree 2 / 4');
A1A3_T2_terms.append(A1A3_graph.boundary_pushforward([-3 * lambdaclass(1, 1, 2), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 4 / 5 of tree 2 / 4');
A1A3_T2_terms.append(A1A3_graph.boundary_pushforward([-1 * fundclass(1, 2), lambdaclass(1, 2, 1), fundclass(1, 1)]));
print('Computed term 5 / 5 of tree 2 / 4');
A1A3_summands.append( 1/1 * R.sum(A1A3_T2_terms));

A1A3_graph = StableGraph([0, 1, 2, 1], [[1, 2, 3], [4], [5], [6]], [(1, 4), (2, 5), (3, 6)]);
A1A3_T3_terms = [];

A1A3_T3_terms.append(A1A3_graph.boundary_pushforward([-3 * fundclass(0, 3), fundclass(1, 1), fundclass(2, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 3 / 4');
A1A3_summands.append( 1/1 * R.sum(A1A3_T3_terms));

A1A3_graph = StableGraph([1, 1, 1, 1], [[1, 2, 3], [4], [5], [6]], [(1, 4), (2, 5), (3, 6)]);
A1A3_T4_terms = [];

A1A3_T4_terms.append(A1A3_graph.boundary_pushforward([1 * fundclass(1, 3), fundclass(1, 1), fundclass(1, 1), fundclass(1, 1)]));
print('Computed term 1 / 1 of tree 4 / 4');
A1A3_summands.append( 1/6 * R.sum(A1A3_T4_terms));

Torelli_pullback = R.sum(A1A3_summands);
