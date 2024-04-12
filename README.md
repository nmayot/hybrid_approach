A git repository that contains scripts in Matlab to performed the analysis presented in Mayot *et al.* (submitted) - Nat. Commun.

## Contact
Nicolas Mayot  
email: n.mayot@uea.ac.uk  

## Folder structure

```text
.
│  
├── /data : the data used are stored here
│   ├── PlankTOM12/: some of the outputs from runs performed with the PlankTOM12.1 model (available to reviewers)
│   ├── zenodo_dataprod/: fCO2-products data, https://zenodo.org/records/10222484
│   ├── zenodo_oceanmodel/: GOBMs data, https://zenodo.org/records/10222484
│   ├── global_co2_merged_2022.txt: xCO2 data provided to ocean modeling teams by GCB in 2022
│   ├── Cflu_RSS_n_RNB.xlsx: results from the hybrid approach (available to reviewers)
│   ├── Global_Carbon_Budget_2023v1.0.xlsx: GCB website, https://globalcarbonbudget.org/carbonbudget2023/
│   ├── EN4_SSS.mat: combined EN.4.2.2 data into a matlab file, https://www.metoffice.gov.uk/hadobs/en4/download-en4-2-2.html
│   ├── RECCAP2_modified.mat: modified reccap2 mask by adding North, Tropics, and South divisions
│   ├── SOCATv2023_tracks_gridded_monthly.nc: SOCAT data, https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0278913
│   ├── sst-OISST1-2.nc: NOAA data, https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html
│   └── surface_air_pressure-ERA5.nc: ERA5 data, https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview
│   
└── /scripts
    ├── calculate_MSE.m   : to read and estimate the values (i.e., CO2 flux and MSE) that will be used by the hybrid approach
    ├── hybrid_approach.m : to perform global, regional, or latitudinal analysis
    ├── fig_XX.m          : to produce a fig XX (1 to 6, and S01 to S06)
    └── table_XX.m        : to calculate the values reported in Table 1 and S02
```

Note that the data are not accessible on GitHub. They need to be downloaded locally from the mentioned sources and placed inside the correct folders. 
