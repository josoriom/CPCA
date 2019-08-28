'use strict';

const { Matrix, EVD, NIPALS } = require('ml-matrix');

// Consensus Principal Component Analysis
class CPCA {
  constructor(x, options = {}) {
    const { condition = 1e-10, blocksize = x.rows } = options;
    // It is the stop condition. When residual matrix norm is equal to 'condition'
    // Block slice (Block-size can be chosen arbitrarily, and this factor will be number of variables by block)
    this.x = x;
    this.x.scale('columns').center('columns');
    let blocks = {};

    let nb = Math.ceil(this.x.columns / blocksize);
    let counter = 0;
    let jump = x.columns;
    for (let i = 0; i < nb; i++) {
      jump -= blocksize;
      let nc = 0;
      if (jump >= 0) {
        nc = blocksize;
      } else {
        nc = blocksize + jump;
      }
      blocks[`block${i}`] = new Matrix(this.x.rows, nc);
      for (let j = 0; j < nc; j++) {
        blocks[`block${i}`].setColumn(j, this.x.getColumnVector(counter));
        counter++;
      }
    }
    // Choosing initial superscore
    let tsup = new EVD(this.x.mmul(this.x.transpose())).V.getColumnVector(0); // The initial tsup is the eigenvector with the largest magnitude eigenvalue of (X.X^t) DOI: 10.1007/s11306-011-0361-9

    tsup = tsup.div(tsup.norm());

    let evbd = [];
    for (let i = 0; i < nb; i++) {
      evbd.push(blocks[`block${i}`].transpose()
        .mmul(blocks[`block${i}`]).trace());
    }
    this.evbd = evbd;

    let Tb = {}; // Object with the blockscores
    let Pb = {}; // Object with blockloadings
    let W = {}; // Object with the blockweights
    let ev = {}; // Object with explained variance by blocks
    let ts = {};
    let test = 1;
    let k = 0;

    do {
      let t = new Matrix(this.x.rows, nb); // Matrix with blockscores by each column
      for (let i = 0; i < nb; i++) {
        t.setColumn(i,
          new NIPALS(blocks[`block${i}`].clone()).t
        );
      }

      ts[`t${k}`] = t;
      this.ts = ts;
      W[`w${k}`] = t.transpose().mmul(tsup);

      tsup = t.mmul(W[`w${k}`]);

      W[`w${k}`] = W[`w${k}`].div(tsup.norm());

      tsup = tsup.div(tsup.norm());

      Tb[`tsup${k}`] = tsup;
      this.Tb = Tb;

      // Deflation
      let p = {};
      let evb = [];
      for (let i = 0; i < nb; i++) {
        p[`pb${i}`] = blocks[`block${i}`].transpose()
          .mmul(tsup);

        blocks[`block${i}`] = blocks[`block${i}`]
          .sub(tsup.mmul(p[`pb${i}`].transpose()));

        evb.push(blocks[`block${i}`].transpose()
          .mmul(blocks[`block${i}`]).trace()
        );
      }

      Pb[`p${k}`] = p;
      ev[`ev${k}`] = evb;
      this.Pb = Pb;
      this.ev = ev;

      let Xr = new Matrix(this.x.rows, this.x.columns); // Residuals matrix of initial data
      counter = 0;
      for (let i = 0; i < nb; i++) {
        for (let j = 0; j < blocks[`block${i}`].columns; j++) {
          Xr.setColumn(counter, blocks[`block${i}`].getColumnVector(j));
          counter++;
        }
      }

      test = Xr.norm();

      tsup = new EVD(Xr.mmul(Xr.transpose())).V.getColumnVector(0);

      tsup = tsup.div(tsup.norm());
      k++;
    } while (test > condition);

    let Tsup = new Matrix(this.x.rows, Object.keys(Tb).length); // Superscores
    let Psup = {}; // Superloadings

    for (let j = 0; j < nb; j++) {
      Psup[`psup${j}`] = new Matrix(Pb[`p${0}`][`pb${j}`].rows, Object.keys(Tb).length);
    }

    for (let i = 0; i < Object.keys(Tb).length; i++) {
      Tsup.setColumn(i, Tb[`tsup${i}`]);
      for (let k = 0; k < nb; k++) {
        Psup[`psup${k}`].setColumn(i, Pb[`p${i}`][`pb${k}`]);
      }
    }

    this.Tsup = Tsup;
    this.Psup = Psup;
  }

  /**
   * Returns a object with the blockscores
   * @returns {[Array]}
   */
  getBlockScores() {
    return this.Tb;
  }

  /**
   * Returns a object with the blockloadings
   * @returns {[Array]}
   */
  getBlockLoadings() {
    return this.Pb;
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
  getSuperLoadings() {
    return this.Psup;
  }

  /**
   * Returns a object with superweights by blocks
   * @returns {[Object]}
   */
  getSuperWeights() {
    let result = {};
    for (let i = 0; i < this.Tsup.columns; i++) {
      let g = [];
      for (let j = 0; j < Object.keys(this.Psup).length; j++) {
        g.push(this.Tsup.getColumnVector(i).transpose()
          .mmul(this.ts[`t${i}`].getColumnVector(j)).get(0, 0)
        );
      }
      result[`w${i}`] = g;
    }
    return result;
  }

  /**
   * Returns a object with explained variance by blocks
   * @returns {[Object]}
   */
  getExplainedVarianceBlocks() {
    let result = {};
    let g = [];
    for (let i = 0; i < Object.keys(this.ev).length; i++) {
      for (let j = 0; j < this.evbd.length; j++) {
        g.push(this.ev[`ev${i}`][j] / this.evbd[j]);
      }
      result[`ev${i}`] = g;
      g = [];
    }
    return result;
  }

  /**
   * Returns a object with amount of variation explained in the data for each superscore
   * @returns {[Object]}
   */
  getLambda() {
    let l = [];
    let s = this.x.mmul(this.x.transpose());
    for (let i = 0; i < this.Tsup.columns; i++) {
      l.push(this.Tsup.getColumnVector(i).transpose().mmul(s)
        .mmul(this.Tsup.getColumnVector(i)).get(0, 0)
      );
    }
    return l;
  }
}

