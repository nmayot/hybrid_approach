A git repository that contains scripts in Matlab to performed the analysis presented in Mayot *et al.* (submitted) - Nat. Commun.

## Contact
Nicolas Mayot  
email: n.mayot@uea.ac.uk  

## Folder structure

```text
.
│  
├── /data : the data used are stored here
│   ├── zenodo_dataprod/  : fCO2-products data, https://zenodo.org/record/7273309
│   ├── global_co2_merged_2022.txt :xCO2 data provided to ocean modeling teams by GCB in 2022
│   ├── Cflu_RSS_n_NoCoast.xlsx : results from the hybrid approach (available to reviewers)
│   ├── Global_Carbon_Budget_2022v1.0.xlsx : GCB website, https://www.globalcarbonproject.org/carbonbudget/22/data.htm
│   ├── EN4_SSS.mat : combined EN.4.2.2 data into a matlab file, https://www.metoffice.gov.uk/hadobs/en4/download-en4-2-2.html
│   ├── RECCAP2_region_masks_all_v20210412.nc : https://reccap2-ocean.github.io/regions/
│   ├── RECCAP2_modified.mat : modified reccap2 mask by adding North, Tropics, South divisions
│   ├── SOCATv2022_tracks_gridded_monthly.nc : SOCAT data, https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0253659
│   ├── sst-OISST1-2.nc : NOAA data, https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
│   └── surface_air_pressure-ERA5.nc : ERA5 data, https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview
│   
└── /scripts
    ├── hybrid_approach.m : to perform global, regional, or latitudinal analysis
    ├── fig_04.m          : to produce a fig 4 (fugacity coeff, open fCO2-products and SOCAT data)
    └── table_01.m        : to calculate the values reported in Table 1
```

Note that the data are not accessible on GitHub. They need to be downloaded locally from the mentioned sources and placed inside the correct folders. 
