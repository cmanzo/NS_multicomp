
# Nested Sampling approach for the multicomponent fitting of data from segmented super-resolution images (NS_multicomp)
We developed a Bayesian approach based on the Nested Sampling algorithm for the multicomponent fitting of data obtained from the segmentation of super-resolution image. The algorithm provides an estimation of the Bayesian evidence to perform model ranking, calculates the model best-fit parameters and their confidence interval. Please refer to the [article](http://dx.doi.org/10.1039/C9CP05616E) below for further information.
A simulated dataset of 300 datapoints with ground truth values *K*=5 and decreasing weights 𝜶=(0.3333,0.2667,0.2000,0.1333,0.0667) is provided for testing. The calibration function used to generate the data is a lognormal distribution with *µ*=3.3491 and 𝜎=0.8462.

## Authors
The R and Matlab code was developed by T. Košuta and C. Manzo in 2019 as part of a research project of [the QuBI lab](https://mon.uvic.cat/qubilab/ "the QuBI lab") carried out at the Faculty of Science and Technology of the Universitat de Vic - Universitat Central de Catalunya (Spain) in collaboration with F. Cella Zanacchi.

## Reference
Please cite the publication below in all your documents and manuscripts that made use of the software included in this repository.
> #### Bayesian analysis of data from segmented super-resolution images for quantifying protein clustering
> ###### Kosuta, T., Cullell-Dalmau, M., Cella Zanacchi, F., Manzo, C.
> ###### *Phys. Chem. Chem. Phys.* 22: 1107-1114 (2020)
> ###### doi: 10.1039/c9cp05616e
> ###### [Link to PhysChemChemPhys](http://dx.doi.org/10.1039/C9CP05616E)
> ###### [Link to ArXiv](https://arxiv.org/abs/1909.13133)

#### *BibTeX*
```
@article{kosuta2020pccp,
  title={Bayesian analysis of data from segmented super-resolution images for quantifying protein clustering},
  author={Ko\v{s}uta, Tina and Cullel-Dalmau, Marta and Cella Zanacchi, Francesca and Manzo, Carlo},
  journal={Phys. Chem. Chem. Phys.},
  volume={22},
  pages={1107--1114},
  year={2020},
  doi={10.1039/c9cp05616e},
  publisher={The Royal Society of Chemistry},
  url  ={http://dx.doi.org/10.1039/C9CP05616E}
}
```
### License
Copyright (c) 2019-2020 Carlo Manzo & Tina Košuta
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
