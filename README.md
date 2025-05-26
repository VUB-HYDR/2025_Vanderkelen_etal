# 2025_Vanderkelen_etal

Python code for lifetime deficit analysis with ISIMIP simulations for "Escalating lifetime water deficit for younger generations". This code is branched off from the Python code of [Grant et al. (2025)](https://github.com/VUB-HYDR/2025_Grant_etal_Nature) and extended toward assessing lifetime water deficits.


## Environment
The python modules used in this repository can be installed using [env_waterscarcity.yml](./env_waterscarcity.yml). This may take up to an hour to compile in the Anaconda prompt. If this doesn't work, you will need to manually collate imported modules in the .py scripts and generate your own env with this list. Try using Python 3.9 for this.

```
conda env create -f env_waterscarcity.yml

```

## Repo structure
The code in this repo is organised as follows: 
1. [preprocessing_isimip](./preprocessing_isimip/) contains the preprocessing code for the ISIMIP data. 
2. [lifetime_deficit](./lifetime_deficit/) contains the code of the lifetime framework and the majority of figures.

See the [overview notebook](./overview.ipynb) for a summary of the code and their purpose. 

## data
Data for our analysis is available [here](). 

## License
This project is licensed under the MIT License. See also the [LICENSE](LICENSE) file.


## References

Grant, L., Vanderkelen, I., Gudmundsson, L., Fischer, E., Seneviratne, S. I., & Thiery, W. (2025). Global emergence of unprecedented lifetime exposure to climate extremes. Nature, 641(8062), 374-379. https://doi.org/10.1038/s41586-025-08907-1


Thiery, W., Lange, S., Rogelj, J., Schleussner, C. F., Gudmundsson, L., Seneviratne, S. I., Andrijevic, M., Frieler, K., Emanuel, K., Geiger, T., Bresch, D. N., Zhao, F., Willner, S. N., Büchner, M., Volkholz, J., Bauer, N., Chang, J., Ciais, P., Dury, M., … Wada, Y. (2021). Intergenerational inequities in exposure to climate extremes. Science, 374(6564), 158–160. https://doi.org/10.1126/science.abi7339
