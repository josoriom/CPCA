import { toBeDeepCloseTo } from 'jest-matcher-deep-close-to';
import { getNumbers } from 'ml-dataset-iris';
import { Matrix } from 'ml-matrix';

import CPCA from '../index';

expect.extend({ toBeDeepCloseTo });

const dataset = getNumbers();
const iris = new Matrix(dataset);

const expectedLoadings = [
  [0.521, 0.269, 0.58, 0.565],
  [0.377, 0.923, 0.024, 0.067],
  [0.72, 0.244, 0.142, 0.634],
  [0.261, 0.124, 0.801, 0.524],
];

describe('iris dataset test method covarianceMatrix', () => {
  it('loadings', () => {
    let cpca = new CPCA(iris, { center: true, scale: true });
    const loadings = cpca
      .getLoadings()
      .eigenVectors.transpose()
      .to2DArray()
      .map((x) => x.map((y) => Math.abs(y)));
    expect(loadings).toBeDeepCloseTo(expectedLoadings, 3);
  });
  it('loadings should be orthogonal', () => {
    let cpca = new CPCA(iris, { center: true, scale: true });
    const m = cpca
      .getLoadings()
      .eigenVectors.transpose()
      .mmul(cpca.getLoadings().eigenVectors)
      .round();
    expect(m.sub(Matrix.eye(4, 4)).sum()).toBe(0);
  });
  it('eigenvalues', () => {
    let cpca = new CPCA(iris, { center: true, scale: true });
    const eigenvalues = cpca.getLoadings().eigenValues;
    expect(eigenvalues).toBeDeepCloseTo(
      [20.853205, 11.67007, 4.676192, 1.756847],
      6,
    );
  });
  it('Relative errors of the aproximation by blocks using 4 components', () => {
    let cpca = new CPCA(iris, {
      center: false,
      scale: false,
      slices: 2,
      components: 4,
    });
    let error = cpca.getDataByBlocks().relativeError;
    expect(error).toBeDeepCloseTo([0.0, 0.0], 2);
  });
});
