import { Matrix, EVD, NIPALS, AbstractMatrix } from 'ml-matrix';

/**
 * Creates new PCA (Principal Component Analysis) from the dataset
 * @param {Matrix} x - dataset
 * @param {Object} [options]
 * @param {number} [options.blocksSlicing] - This option indicates how the partition of the blocks will be done. Could be a 'number' in the case of data blocks with equal number of blocks or 'Array' in the case of data blocks with different number of columns. In the 'Array' case each element of the array indicates one block and the number indicates the number of columns in the corresponding block (default: number of variables in the data set).
 * @param {boolean} [options.componentsNumber] - This option indicates the number of components that will work with (default: number of variables in the data set).
 * @param {boolean} [options.center] - boolean indicating if data should be centered (default: 'true').
 * @param {boolean} [options.scale] - boolean indicating if data should be scaled (default: 'false').
 * @param {number} [options.method] - It refers to the method to calculate the initial super score. This method can be 'EVD' or 'ones' (default: 'ones').
 * */

export default class CPCA {
  public center?: boolean;
  public scale?: boolean;
  public slices: number[];
  public components?: number;
  public method?: string;
  public data: Matrix;
  public loadings: Matrix[][];
  public weigths: Matrix[];
  public superScores: Matrix;
  constructor(x: Matrix, options: Options = {}) {
    x = x.clone();
    this.data = x;
    const {
      center = true,
      scale = false,
      method = 'ones',
      slices = [x.columns],
      components = x.columns,
    } = options;

    if (center) {
      x.center('column');
    }

    if (scale) {
      x.scale('column');
    }

    let P: Matrix[][] = [];
    let W: Matrix[] = [];

    let blocks: Matrix[] = [];
    this.components = components;
    if (Array.isArray(slices)) {
      blocks = [];
      this.slices = slices;
      let counter = 0;
      for (let slice of slices) {
        let h = new Matrix(x.rows, slice);
        for (let j = 0; j < slice; j++) {
          h.setColumn(j, x.getColumnVector(counter));
          counter++;
        }
        blocks.push(h);
      }
    } else if (typeof slices === 'number') {
      let variablesByBlock: number[] = [];
      blocks = [];
      let counter = 0;
      let jump = x.columns;
      for (let i = 0; i < Math.ceil(x.columns / slices); i++) {
        jump -= slices;
        let nc = 0;
        if (jump >= 0) {
          nc = slices;
        } else {
          nc = slices + jump;
        }
        variablesByBlock.push(nc);
        let h = new Matrix(x.rows, nc);
        for (let j = 0; j < nc; j++) {
          h.setColumn(j, x.getColumnVector(counter));
          counter++;
        }
        blocks.push(h);
      }
      this.slices = variablesByBlock;
    } else {
      this.slices = slices;
    }

    // Choosing initial superscore
    let tsup:
      | ArrayLike<number>
      | (AbstractMatrix | ArrayLike<ArrayLike<number>>);
    if (method === 'EVD') {
      tsup = new EVD(x.mmul(x.transpose())).eigenvectorMatrix.getColumnVector(
        0,
      );
    } else {
      tsup = new Matrix(x.rows, 1).setColumn(0, new Array(x.rows).fill(1));
    }
    tsup = tsup.div(tsup.norm('frobenius'));
    let Tsup = new Matrix(x.rows, components);
    let k = 0;
    do {
      let T = new Matrix(x.rows, blocks.length);
      let w: Matrix;
      for (let i = 0; i < blocks.length; i++) {
        let a = blocks[i].clone();
        let h = new NIPALS(a);
        T.setColumn(i, h.t.div(Math.sqrt(blocks[i].columns)));
      }
      w = T.transpose().mmul(tsup);
      tsup = T.mmul(w);
      W.push(w.div(w.norm('frobenius')));
      tsup = tsup.div(tsup.norm('frobenius'));
      Tsup.setColumn(k, tsup);
      // Deflation
      let p: Matrix[] = [];
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
        for (let block of blocks) {
          for (let j = 0; j < block.columns; j++) {
            Xr.setColumn(counter, block.getColumnVector(j));
            counter++;
          }
        }
        tsup = new EVD(
          Xr.mmul(Xr.transpose()),
        ).eigenvectorMatrix.getColumnVector(0);
      } else if (method === 'ones') {
        tsup = new Matrix(x.rows, 1).setColumn(0, new Array(x.rows).fill(1));
      }
      tsup = tsup.div(tsup.norm('frobenius'));
      k++;
    } while (k < components);
    this.weigths = W;
    this.superScores = Tsup;
    this.loadings = P;
  }

  getSuperScores() {
    return this.superScores;
  }

  /**
   * Returns a object with superloadings by blocks
   */
  getLoadings() {
    let loadingBlocks: Matrix[] = []; // Building Loadings by blocks
    for (let i = 0; i < this.loadings[0].length; i++) {
      let h = new Matrix(this.loadings[0][i].rows, this.components as number);
      for (let j = 0; j < (this.components as number); j++) {
        h.setColumn(j, this.loadings[j][i]);
      }
      loadingBlocks.push(h);
    }
    // Building Superloadings
    let Psup = new Matrix(this.data.columns, this.components as number);
    let counter = 0;
    for (let i = 0; i < this.loadings[0].length; i++) {
      if (i === 0) {
        Psup.setSubMatrix(loadingBlocks[i], counter, 0);
      } else if (this.slices[i] === 1) {
        Psup.setSubMatrix(loadingBlocks[i], counter, 0);
      }
      counter += loadingBlocks[i].rows;
    }

    let eigenValues: number[] = [];
    let eigenVectors = new Matrix(Psup.rows, Psup.columns);
    for (let i = 0; i < Psup.columns; i++) {
      let z = Psup.getColumnVector(i).norm('frobenius');
      eigenVectors.setColumn(i, Psup.getColumnVector(i).div(z));
      eigenValues.push(z);
    }
    let result = { loadingBlocks, Psup, eigenValues, eigenVectors };
    return result;
  }

  /**
   * Returns a object with superweights by blocks
   */
  getSuperWeights() {
    return this.weigths;
  }

  /**
   * Returns a object with amount of variation explained in the data for each superscore
   */
  getDataByBlocks() {
    let predictedBlocks: Matrix[] = [];
    let dataBlocks: Matrix[] = [];
    let counter = 0;
    let relativeError: number[] = [];
    for (let j = 0; j < this.loadings[0].length; j++) {
      let h: Matrix;
      let g = new Matrix(this.data.rows, this.slices[j]);
      for (let i = 0; i < this.loadings.length; i++) {
        h = this.superScores
          .getColumnVector(i)
          .mmul(this.loadings[i][j].clone().transpose());
        g.add(h);
      }
      predictedBlocks.push(g);
      let data = new Matrix(this.data.rows, g.columns);
      for (let i = 0; i < g.columns; i++) {
        data.setColumn(i, this.data.getColumnVector(counter));
        counter++;
      }
      dataBlocks.push(data);
      let factor =
        dataBlocks[j].norm('frobenius') / predictedBlocks[j].norm('frobenius');
      let s = dataBlocks[j].clone().sub(predictedBlocks[j].clone().mul(factor));
      let o = s.norm('frobenius') / dataBlocks[j].norm('frobenius');
      relativeError.push(o);
    }
    let result = { relativeError, predictedBlocks };
    return result;
  }
}

/**
 * CPCA options
 */
interface Options {
  center?: boolean;
  scale?: boolean;
  slices?: number;
  components?: number;
  method?: string;
}
