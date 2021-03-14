![image](https://github.com/patternizer/glosat-20CRv3-t2m/blob/master/t2m_1815_01.png)

# glosat-20CRv3-t2m

SLURM shell templates and python script generator code to download and process monthly 80-member ensemble 2m-temperature data from 20CRv3 reanalysis and compare with co-located GloSAT.p0x station land surface air temperature data as part of ongoing work for the [GloSAT](https://www.glosat.org) project: www.glosat.org. 

## Contents

* `SLURM/` - shell scripts and python script generator codes for the initial processing chain 
* `extract_timeseries_runnable_absolutes.py` - python script to extract 20CRv3 2m absolute temperatures and calculate monthly ensemble averages and percentiles for a specific gridcell (lat,lon)
* `extract_timeseries_runnable_anomalies.py` - python script to extract 20CRv3 absolute temperatures and calculate monthly ensemble averages and percentiles for a specific gridcell (lat,lon) relative to the 1961-1990 baseline mean
* `extract_timeseries_and_compare_absolutes.py` - python script to extract the 20CRv3 monthly mean absolute temperature timeseries co-located with GloSAT.p03 station data
* `extract_timeseries_and_compare_anomalies.py` - python script to extract the 20CRv3 monthly mean absolute temperature timeseries co-located with GloSAT.p03 station data relative to the 1961-1990 baseline mean
* `plot_ensemble_mean_monthly_anomalies_contactsheet.py` - python script to plot yearly contact sheets of 12-monthly global ensemble mean maps in Orthographic projection (centred on the Greenwich meridian)
* `plot_ensemble_mean_monthly_anomalies_months1-12.py` - python script to plot monthly ensemble mean maps in Robinson projection (centred on the Greenwich meridian)

The first step is to clone the latest glosat-20CRv3-t2m code into your /home/users/ folder on JASMIN and step into the installed Github directory: 

    $ git clone https://github.com/patternizer/glosat-20CRv3-t2m.git
    $ cd glosat-20CRv3-t2m

### Using Standard Python

The code is designed to run on CentOS-7 systems with the SLURM batch process scheduler and the JasPy module.
Processing chain code for new implementations is in /SLURM.

* Example comparison runs:

First lon on and navigate to the codebase in the GloSAT GWS at /development/data/raw/UEA/20CRv3/t2m/monthly/ and invoke:

    $ module load jaspy 

Extract 20CRv3 2m temperature ensemble averages and percentiles from a specific gridcell (lat,lon) as follows:

    $ python extract_timeseries_runnable_absolutes.py -- lat lon
    $ python extract_timeseries_runnable_anomalies.py -- lat lon

The '--' is required to overide options in the case of negative valued lat or lon coordinates.
For example:

    $ python extract_timeseries_runnable_absolutes.py 52.5 13.3 (OR)
    $ python extract_timeseries_runnable_absolutes.py -- 52.5 13.3
    $ python extract_timeseries_runnable_anomalies.py -- -35.9 150.2

These calls generate a plot of the 30-year rolling average ensemble mean and median as well as the 5-95,10-90 and 25-75 percentile ranges for the extracted ensemble absolute temperature or anomalies (from 1961-1990) for the gridcell containing (lat,lon). The Haversine distance in km to the centre of the gridcell is reported. 
For example, for extracted anomalies at (55.5°N,-60.2°E):

![image](https://github.com/patternizer/glosat-20CRv3-t2m/blob/master/lon_299.8_lat_55.5.png)

Both the figure and CSV files containing the monthly extracted statistics timeseries will be output to the folder EXTRACT_RUN/.

To compare with co-located GloSAT.p03 station data, extract 20CRv3 mean absolute temperature or anomalies (from 1961-1990) by adding either a unique search term from a station name or its station code as follows:

    $ python extract_timeseries_and_compare_absolutes.py oxford
    $ python extract_timeseries_and_compare_absolutes.py 038900
    $ python extract_timeseries_and_compare_anomalies.py uppsala
    $ python extract_timeseries_and_compare_anomalies.py 024581

This will generate a plot of the yearly averaged extracted 20CRv3 t2m timeseries together with the available station timeseries.
For example, for Jan Mayen station anomalies at (70.9°N,-8.7°E):

![image](https://github.com/patternizer/glosat-20CRv3-t2m/blob/master/010010_jan_mayen_anomaly_yearly.png)

* Example summary plots:

    $ python plot_ensemble_mean_monthly_anomalies_contactsheet.py

This will generate a yearly contact sheet of 12 monthly global ensemble mean maps in Orthographic projection (centred on the Greenwich meridian).
For example, for 2015:

![image](https://github.com/patternizer/glosat-20CRv3-t2m/blob/master/t2m_2015.png)

    $ python plot_ensemble_mean_monthly_anomalies_months1-12.py

This will generate a monthly ensemble mean map in Robinson projection (centred on the Greenwich meridian).
For example, for January 2015:

![image](https://github.com/patternizer/glosat-20CRv3-t2m/blob/master/t2m_2015_01.png)

* Animating summary plots:

An animated GIF can be generated with:

    $ convert -delay 100 -loop 0 t2m_*.png t2m.gif

and converted to MP4 with:

    $ ffmpeg -i t2m.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" t2m.mp4
    
## License

The code is distributed under terms and conditions of the [Open Government License](http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

## Contact information

* [Michael Taylor](michael.a.taylor@uea.ac.uk)

