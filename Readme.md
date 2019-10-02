# ml-cpca

  [![NPM version][npm-image]][npm-url]
  [![build status][travis-image]][travis-url]
  [![David deps][david-image]][david-url]
  [![npm download][download-image]][download-url]
  
Consensus Principal Component Analysis

# CPCA
<p align="center">
</p>
<p align="center">
  <img alt="CPCA src="images/cpca.png">
</p>

## Installation

`$ npm install ml-cpca`

## [API Documentation](https://cheminfo.github.io/cpca/)

## Example

```js
const cpca = require('ml-cpca');
```
**Equations**
* Data set partition where nb is the number of blocks and N is the number of components.

![equation](https://tex.cheminfo.org/?tex=X%20=%20\left%20[%20X_{1},%20...,%20X_{n_{b}}%20\right%20])

* Super score (Tsup) is a Matrix with block scores (tsup,j) arranged by columns.

![equation](https://tex.cheminfo.org/?tex=T_{sup}%20=%20[t_{sup,1},%20...,%20t_{sup,N}])

* Super loadings (Psup) is an object with block loadings (psup,j).

![equation](https://tex.cheminfo.org/?tex=P_{sup}%20=%20[p_{sup,1},%20...,%20p_{sup,N}])


* Each block can be rebuilded by the outer product of its respective block score and block loading.

![equation](https://tex.cheminfo.org/?tex=X_{k}%20=%20\sum_{j=1}^{N}t_{sup,k}^{j}\otimes%20p_{sup,k}^{j}%20%20;%20for\\%20%20%20k%20=%201,%20...,%20n_{b})

**Arguments**

* `x`: Matrix containing the inputs.

**Options**

* `center`: boolean indicating if data should be centered (default: 'true').
* `scale`: boolean indicating if data should be scaled (default: 'false').
* `method`: It refers to the method to calculate the initial super score. This method can be 'EVD' or 'ones' (default: 'ones').
* `blocksSlicing`: This option indicates how the partition of the blocks will be done. Could be a 'number' in the case of data blocks with equal number of blocks or 'Array' in the case of data blocks with different number of columns. In the 'Array' case each element of the array indicates one block and the number indicates the number of columns in the corresponding block (default: number of variables in the data set).
* `componentsNumber`: This option indicates the number of components that will work with (default: number of variables in the data set).

## License

[MIT](./LICENSE)

[npm-image]: https://img.shields.io/npm/v/cpca.svg?style=flat-square
[npm-url]: https://www.npmjs.com/package/cpca
[travis-image]: https://img.shields.io/travis/cheminfo/cpca/master.svg?style=flat-square
[travis-url]: https://travis-ci.org/cheminfo/cpca
[david-image]: https://img.shields.io/david/cheminfo/cpca.svg?style=flat-square
[david-url]: https://david-dm.org/cheminfo/cpca
[download-image]: https://img.shields.io/npm/dm/cpca.svg?style=flat-square
[download-url]: https://www.npmjs.com/package/cpca

## References

[1] Westerhuis, J. A., Kourti, T., & MacGregor, J. F. (1998). Analysis of multiblock and hierarchical PCA and PLS models. Journal of Chemometrics: A Journal of the Chemometrics Society, 12(5), 301-321. DOI: https://doi.org/10.1002/(SICI)1099-128X(199809/10)12:5<301::AID-CEM515>3.0.CO;2-S

[2] Smilde, A. K., Westerhuis, J. A., & de Jong, S. (2003). A framework for sequential multiblock component methods. Journal of Chemometrics: A Journal of the Chemometrics Society, 17(6), 323-337. DOI: https://doi.org/10.1002/cem.811

[3] Xu, Y., & Goodacre, R. (2012). Multiblock principal component analysis: an efficient tool for analyzing metabolomics data which contain two influential factors. Metabolomics, 8(1), 37-51. DOI: https://doi.org/10.1007/s11306-011-0361-9
