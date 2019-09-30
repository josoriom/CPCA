'use strict';

const { Matrix, EVD, NIPALS } = require('ml-matrix');

/**
 * Creates new PCA (Principal Component Analysis) from the dataset
 * @param {Matrix} x - dataset
 * @param {Object} [options]
 * @param {numeric} [options.blocksSlicing] - This option indicates how the partition of the blocks will be done. Could be a 'number' in the case of data blocks with equal number of blocks or 'Array' in the case of data blocks with different number of columns. In the 'Array' case each element of the array indicates one block and the number indicates the number of columns in the corresponding block (default: number of variables in the data set).
 * @param {boolean} [options.componentsNumber] - This option indicates the number of components that will work with (default: number of variables in the data set).
 * @param {boolean} [options.center] - boolean indicating if data should be centered (default: 'true').
 * @param {boolean} [options.scale] - boolean indicating if data should be scaled (default: 'false').
 * @param {numeric} [options.method] - It refers to the method to calculate the initial super score. This method can be 'EVD' or 'ones' (default: 'ones').
 * */

class CPCA {
  constructor(x, options = {}) {
    this.x = x.clone();
    if (x === true) {
      this.center = options.center;
      this.scale = options.scale;
      this.blocksSlicing = options.blocksSlicing;
      this.componentsNumber = componentsNumber;
      this.method = options.method;
      return;
    }
    const {
      center = true,
      scale = false,
      method = 'ones',
      blocksSlicing = x.columns,
      componentsNumber = x.columns
    } = options;
    if (center === true) {
      this.x.center('column');
    } if (scale === true) {
      this.x.scale('column');
    }

    let blocks, numberOfblocks, counter;
    this.componentsNumber = componentsNumber;
    if (Array.isArray(blocksSlicing)) {
      this.blocksSlicing = blocksSlicing;
      blocks = {};
      numberOfblocks = blocksSlicing.length;
      counter = 0;
      for (let i = 0; i < blocksSlicing.length; i++) {
        blocks[`block${i}`] = new Matrix(this.x.rows, blocksSlicing[i]);
        for (let j = 0; j < blocksSlicing[i]; j++) {
          blocks[`block${i}`].setColumn(j, this.x.getColumnVector(counter));
          counter++;
        }
      }
    } if (typeof (blocksSlicing) === 'number') {
      let variablesByBlock = [];
      blocks = {};
      numberOfblocks = Math.ceil(this.x.columns / blocksSlicing);
      counter = 0;
      let jump;
      jump = x.columns;
      for (let i = 0; i < numberOfblocks; i++) {
        jump -= blocksSlicing;
        let nc = 0;
        if (jump >= 0) {
          nc = blocksSlicing;
        } else {
          nc = blocksSlicing + jump;
        }
        variablesByBlock.push(nc);
        blocks[`block${i}`] = new Matrix(this.x.rows, nc);
        for (let j = 0; j < nc; j++) {
          blocks[`block${i}`].setColumn(j, this.x.getColumnVector(counter));
          counter++;
        }
      }
      this.blocksSlicing = variablesByBlock;
    }

    // Choosing initial superscore
    let tsup;
    if (method === 'EVD') {
      tsup = new EVD(this.x.mmul(this.x.transpose())).V.getColumnVector(0);
    } if (method === 'ones') {
      tsup = new Matrix(this.x.rows, 1).setColumn(0, new Array(this.x.rows).fill(1));
    }
    tsup = tsup.div(tsup.norm());
    let Tsup = new Matrix(this.x.rows, componentsNumber);
    let P = {};
    let W = {};
    let k = 0;
    do {
      let T = new Matrix(this.x.rows, numberOfblocks);
      let w;
      for (let i = 0; i < numberOfblocks; i++) {
        let g = blocks[`block${i}`].clone();
        let h = new NIPALS(g);
        T.setColumn(i, h.t.div(Math.sqrt(blocks[`block${i}`].columns)));
      }
      w = T.transpose().mmul(tsup);
      tsup = T.mmul(w);
      W[`w${k}`] = w.div(w.norm());
      tsup = tsup.div(tsup.norm());
      Tsup.setColumn(k, tsup);
      // Deflation
      let p = {};
      for (let i = 0; i < numberOfblocks; i++) {
        p[`block${i}`] = blocks[`block${i}`].transpose()
          .mmul(tsup);
        blocks[`block${i}`] = blocks[`block${i}`]
          .sub(tsup.mmul(p[`block${i}`].transpose()));
      }
      P[`component${k}`] = p;
      // Updating tsup
      if (method === 'EVD') {
        let Xr = new Matrix(this.x.rows, this.x.columns);
        counter = 0;
        for (let i = 0; i < numberOfblocks; i++) {
          for (let j = 0; j < blocks[`block${i}`].columns; j++) {
            Xr.setColumn(counter, blocks[`block${i}`].getColumnVector(j));
            counter++;
          }
        }
        tsup = new EVD(Xr.mmul(Xr.transpose())).V.getColumnVector(0);
      } if (method === 'ones') {
        tsup = new Matrix(x.rows, 1)
          .setColumn(0, new Array(x.rows).fill(1));
      }
      tsup = tsup.div(tsup.norm());
      k++;
    } while (k < componentsNumber);
    this.W = W;
    this.Tsup = Tsup;
    this.P = P;
  }

  /**
   * Returns a matrix with superscores by columns
   * @returns {[Array]}
   */
  getSuperScores() {
    return this.Tsup;
  }

  /**
   * Returns a object with superloadings by blocks
   * @returns {[Object]}
   */
  getLoadings() {
    let loadingBlocks = {}; // Building Loadings by blocks
    for (let i = 0; i < Object.keys(this.P.component0).length; i++) {
      loadingBlocks[`loadingBlock_${i}`] = new Matrix(this.P[`component${0}`][`block${i}`].rows, this.componentsNumber);
      for (let j = 0; j < this.componentsNumber; j++) {
        loadingBlocks[`loadingBlock_${i}`].setColumn(j, this.P[`component${j}`][`block${i}`]);
      }
    }
    // Building Superloadings
    let Psup = new Matrix(this.x.columns, this.componentsNumber);
    for (let i = 0; i < Object.keys(this.P.component0).length; i++) {
      if (i === 0) {
        Psup.setSubMatrix(loadingBlocks[`loadingBlock_${i}`], 0, 0);
      } else {
        if (this.blocksSlicing === 1) {
          Psup.setSubMatrix(loadingBlocks[`loadingBlock_${i}`], i, 0);
        } else {
          Psup.setSubMatrix(loadingBlocks[`loadingBlock_${i}`], loadingBlocks[`loadingBlock_${i - 1}`].rows, 0);
        }
      }
    }

    let eigenValues = []; let eigenVectors = new Matrix(Psup.rows, Psup.columns);
    for (let i = 0; i < Psup.columns; i++) {
      let z = Psup.getColumnVector(i).norm();
      eigenVectors.setColumn(i, Psup.getColumnVector(i).div(z));
      eigenValues.push(z);
    }
    let result = { loadingBlocks, Psup, eigenValues, eigenVectors };
    return result;
  }

  /**
   * Returns a object with superweights by blocks
   * @returns {[Object]}
   */
  getSuperWeights() {
    return this.W;
  }

  /**
   * Returns a object with amount of variation explained in the data for each superscore
   * @returns {[Object]}
   */
  getDataByBlocks() {
    let predictedBlocks = {}; let dataBlocks = {}; let counter = 0; let relativeError = [];
    for (let j = 0; j < Object.keys(this.P.component0).length; j++) {
      let h, g;
      g = new Matrix(this.x.rows, this.blocksSlicing[j]);
      for (let i = 0; i < Object.keys(this.P).length; i++) {
        h = this.Tsup.getColumnVector(i)
          .mmul(this.P[`component${i}`][`block${j}`].clone().transpose());
        g.add(h);
      }
      predictedBlocks[`Block_${j}`] = g;
      let data = new Matrix(this.x.rows, g.columns);
      for (let i = 0; i < g.columns; i++) {
        data.setColumn(i, this.x.getColumnVector(counter));
        counter++;
      }
      dataBlocks[`Block_${j}`] = data;
      let factor = dataBlocks[`Block_${j}`].norm() / predictedBlocks[`Block_${j}`].norm();
      let s = dataBlocks[`Block_${j}`].clone().sub(predictedBlocks[`Block_${j}`].clone().mul(factor));
      let o = s.norm() / dataBlocks[`Block_${j}`].norm();
      relativeError.push(o);
    }
    let result = { relativeError, predictedBlocks };
    return result;
  }
}
module.exports = CPCA;

