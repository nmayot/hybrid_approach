A git repository that contains scripts in Matlab to performed the analysis presented in Mayot *et al.* (in revision) - Nat. Commun.

## Contact
Nicolas Mayot  
email: nicolas.mayot@imev-mer.fr  

## Folder structure

```text
.
│  
├── /data : the data used are stored here
│   ├── PlankTOM12/: outputs from PlankTOM12.1 (https://osf.io/2kzps/?view_only=6ad809f1887342a0a19907e40a33e7cf)
│   ├── zenodo_dataprod/: fCO2-products data, https://zenodo.org/records/10222484
│   ├── zenodo_oceanmodel/: GOBMs data, https://zenodo.org/records/10222484
│   ├── global_co2_merged_2022.txt: xCO2 data provided by GCB in 2022 (https://osf.io/2kzps/?view_only=6ad809f1887342a0a19907e40a33e7cf)
│   ├── Cflu_RSS_n_RNB.xlsx: results from the hybrid approach (https://osf.io/2kzps/?view_only=6ad809f1887342a0a19907e40a33e7cf)
│   ├── Global_Carbon_Budget_2023v1.0.xlsx: GCB website (https://osf.io/2kzps/?view_only=6ad809f1887342a0a19907e40a33e7cf)
│   ├── EN4_SSS.mat: combined EN.4.2.2 data (https://www.metoffice.gov.uk/hadobs/en4/download-en4-2-2.html) into a matlab file (https://osf.io/2kzps/?view_only=6ad809f1887342a0a19907e40a33e7cf)
│   ├── RECCAP2_modified.mat: modified reccap2 mask (https://osf.io/2kzps/?view_only=6ad809f1887342a0a19907e40a33e7cf)
│   ├── SOCATv2023_tracks_gridded_monthly.nc: SOCAT data (https://socat.info/index.php/data-access/)
│   ├── sst-OISST1-2.nc: NOAA data (https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html)
│   ├── dataprod_GCB_19-23.xlsx: data for figure 6 (https://osf.io/2kzps/?view_only=6ad809f1887342a0a19907e40a33e7cf)
│   ├── GCB2023_metrics.xlsx: data for supplementary figures (https://osf.io/2kzps/?view_only=6ad809f1887342a0a19907e40a33e7cf)
│   └── surface_air_pressure-ERA5.nc: ERA5 data, (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview)
│   
└── /scripts
    ├── calculate_MSE.m   : to read and estimate the values (i.e., CO2 flux and MSE) that will be used by the hybrid approach
    ├── hybrid_approach.m : to perform global, regional, or latitudinal analysis
    ├── calculate_dfCO2.m : values that are used in figure 4
    ├── fig_XX.m          : to produce a fig XX (1 to 6, and S01 to S06)
    └── table_XX.m        : to calculate the values reported in Table 1 and S02
```

Note that the data are not accessible on GitHub. They need to be downloaded locally from the mentioned sources and placed inside the correct folders. 
