# size-splash

Collection of scripts used to measure the mass-size relation of galaxy clusters from the Hydrangea simulations, see Contigiani+ 2020 ([arXiv](https://arxiv.org/abs/2012.01336)).

## Files

* load.py and load_galaxies.py and load_mass_profiles.py : loops designed to load Hydrangea hd5 files and write relevant data into numpy bin files. 

* cluster.py : handy class representing the individual Hydrangea clusters and their properties. 

* profiles.py : mass/density profiles used in the paper. 

## Jupyter Notebooks

* analysis.ipynb : notebook where the loaded data is further analyzed.

* plot.ipynb : interactive notebook used to generate the figures found in the paper.


## References

* What is splashback? Benedikt Diemer and Andrey V. Kravtsov, 2014, APJ ([arXiv](https://arxiv.org/abs/1401.1216))

* Splashback measurement : Contigiani+ 2018 ([arXiv](http://arxiv.org/abs/1702.01722), [GitHub](https://github.com/contigiani/splash))

* Everything here is expected to work on the data presented in the original Hydrangea paper ([arXiv](https://arxiv.org/abs/1703.10610))