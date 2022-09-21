import { Matrix, EVD, NIPALS } from 'ml-matrix';

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

export default class CPCA {
  constructor(x, options = {}) {
    x = x.clone();
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
      componentsNumber = x.columns,
    } = options;
    if (center === true) {
      x.center('column');
    }

    if (scale === true) {
      x.scale('column');
    }

    let blocks;
    this.componentsNumber = componentsNumber;
    if (Array.isArray(blocksSlicing)) {
      blocks = [];
      this.blocksSlicing = blocksSlicing;
      let counter = 0;
      for (let i = 0; i < blocksSlicing.length; i++) {
        let h = new Matrix(x.rows, blocksSlicing[i]);
        for (let j = 0; j < blocksSlicing[i]; j++) {
          h.setColumn(j, x.getColumnVector(counter));
          counter++;
        }
        blocks.push(h);
      }
    } else if (typeof blocksSlicing === 'number') {
      let variablesByBlock = [];
      blocks = [];
      let counter = 0;
      let jump = x.columns;
      for (let i = 0; i < Math.ceil(x.columns / blocksSlicing); i++) {
        jump -= blocksSlicing;
        let nc = 0;
        if (jump >= 0) {
          nc = blocksSlicing;
        } else {
          nc = blocksSlicing + jump;
        }
        variablesByBlock.push(nc);
        let h = new Matrix(x.rows, nc);
        for (let j = 0; j < nc; j++) {
          h.setColumn(j, x.getColumnVector(counter));
          counter++;
        }
        blocks.push(h);
      }
      this.blocksSlicing = variablesByBlock;
    }

    // Choosing initial superscore
    let tsup;
    if (method === 'EVD') {
      tsup = new EVD(x.mmul(x.transpose())).V.getColumnVector(0);
    } else if (method === 'ones') {
      tsup = new Matrix(x.rows, 1).setColumn(0, new Array(x.rows).fill(1));
    }
    tsup = tsup.div(tsup.norm());
    let Tsup = new Matrix(x.rows, componentsNumber);
    let P = [];
    let W = [];
    let k = 0;
    do {
      let T = new Matrix(x.rows, blocks.length);
      let w;
      for (let i = 0; i < blocks.length; i++) {
        let a = blocks[i].clone();
        let h = new NIPALS(a);
        T.setColumn(i, h.t.div(Math.sqrt(blocks[i].columns)));
      }
      w = T.transpose().mmul(tsup);
      tsup = T.mmul(w);
      W.push(w.div(w.norm()));
      tsup = tsup.div(tsup.norm());
      Tsup.setColumn(k, tsup);
      // Deflation
      let p = [];
      for (let i = 0; i < blocks.length; i++) {
        let g = blocks[i].transpose().mmul(tsup);
        blocks[i] = blocks[i].sub(tsup.mmul(g.transpose()));
        p.push(g);
      }
      P.push(p);
      // Updating tsup
      if (method === 'EVD') {
        let Xr = new Matrix(x.rows, x.columns);
        let counter = 0;
        for (let i = 0; i < blocks.length; i++) {
          for (let j = 0; j < blocks[i].columns; j++) {
            Xr.setColumn(counter, blocks[i].getColumnVector(j));
            counter++;
          }
        }
        tsup = new EVD(Xr.mmul(Xr.transpose())).V.getColumnVector(0);
      } else if (method === 'ones') {
        tsup = new Matrix(x.rows, 1).setColumn(0, new Array(x.rows).fill(1));
      }
      tsup = tsup.div(tsup.norm());
      k++;
    } while (k < componentsNumber);
    this.x = x;
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
    let loadingBlocks = []; // Building Loadings by blocks
    for (let i = 0; i < this.P[0].length; i++) {
      let h = new Matrix(this.P[0][i].rows, this.componentsNumber);
      for (let j = 0; j < this.componentsNumber; j++) {
        h.setColumn(j, this.P[j][i]);
      }
      loadingBlocks.push(h);
    }
    // Building Superloadings
    let Psup = new Matrix(this.x.columns, this.componentsNumber);
    let counter = 0;
    for (let i = 0; i < this.P[0].length; i++) {
      if (i === 0) {
        Psup.setSubMatrix(loadingBlocks[i], counter, 0);
      } else if (this.blocksSlicing[i] === 1) {
        Psup.setSubMatrix(loadingBlocks[i], counter, 0);
      }
      counter += loadingBlocks[i].rows;
    }

    let eigenValues = [];
    let eigenVectors = new Matrix(Psup.rows, Psup.columns);
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
    let predictedBlocks = [];
    let dataBlocks = [];
    let counter = 0;
    let relativeError = [];
    for (let j = 0; j < this.P[0].length; j++) {
      let h, g;
      g = new Matrix(this.x.rows, this.blocksSlicing[j]);
      for (let i = 0; i < this.P.length; i++) {
        h = this.Tsup.getColumnVector(i).mmul(this.P[i][j].clone().transpose());
        g.add(h);
      }
      predictedBlocks.push(g);
      let data = new Matrix(this.x.rows, g.columns);
      for (let i = 0; i < g.columns; i++) {
        data.setColumn(i, this.x.getColumnVector(counter));
        counter++;
      }
      dataBlocks.push(data);
      let factor = dataBlocks[j].norm() / predictedBlocks[j].norm();
      let s = dataBlocks[j].clone().sub(predictedBlocks[j].clone().mul(factor));
      let o = s.norm() / dataBlocks[j].norm();
      relativeError.push(o);
    }
    let result = { relativeError, predictedBlocks };
    return result;
  }
}
