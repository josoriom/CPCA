'use strict';

const { Matrix } = require('ml-matrix');
const { getNumbers } = require('ml-dataset-iris');
const { toBeDeepCloseTo } = require('jest-matcher-deep-close-to');

const CPCA = require('../index');

expect.extend({ toBeDeepCloseTo });

const iris = new Matrix(getNumbers());

const expectedLoadings = [
  [0.521, 0.269, 0.58, 0.565],
  [0.377, 0.923, 0.024, 0.067],
  [0.72, 0.244, 0.142, 0.634],
  [0.261, 0.124, 0.801, 0.524],
];

describe('iris dataset test method covarianceMatrix', function () {
  var cpca = new CPCA(iris, { center: true, scale: true });
  it('loadings', function () {
    var loadings = cpca
      .getLoadings()
      .eigenVectors
      .transpose()
      .to2DArray()
      .map((x) => x.map((y) => Math.abs(y)));
    expect(loadings).toBeDeepCloseTo(expectedLoadings, 3);
  });
  it('loadings should be orthogonal', function () {
    let m = cpca
      .getLoadings()
      .eigenVectors
      .transpose()
      .mmul(cpca.getLoadings().eigenVectors)
      .round();
    expect(m.sub(Matrix.eye(4, 4)).sum()).toStrictEqual(0);
  });

  it('eigenvalues', function () {
    let eigenvalues = cpca.getLoadings().eigenValues;
    expect(eigenvalues).toBeDeepCloseTo(
      [20.853205, 11.67007, 4.676192, 1.756847],
      6,
    );
  });
});

describe('First block of iris dataset from blockloadings and blockscores', function () {
  it('loadingsh', function () {
    var data = new CPCA(iris, { center: false, scale: false, blocksSlicing: 2 });
    var blocks = data
      .getDataByBlocks()
      .predictedBlocks
      .Block_0
      .getRowVector(0)
      .to2DArray()[0];
    expect(blocks).toBeDeepCloseTo([5.1, 3.5], 2);
  });
  it('Second block of iris dataset from blockloadings and blockscores', function () {
    var data = new CPCA(iris, { center: false, scale: false, blocksSlicing: 2 });
    var blocks = data
      .getDataByBlocks()
      .predictedBlocks
      .Block_1
      .getRowVector(0)
      .to2DArray()[0];
    expect(blocks).toBeDeepCloseTo([1.4, 0.2], 2);
  });
  it('Relative errors of the aproximation by blocks using 4 components', function () {
    var data = new CPCA(iris, { center: false, scale: false, blocksSlicing: 2, componentsNumber: 4 });
    var error = data
      .getDataByBlocks()
      .relativeError;
    expect(error).toBeDeepCloseTo([0.0, 0.0], 2);
  });
});
