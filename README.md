# jasmin-20CRv3

SLURM shell templates including python script generator code to download monthly 80-member ensemble variables (2m Temperature used as an example) from 20CRv3 reanalysis as part of ongoing work for the [GloSAT](https://www.glosat.org) project: www.glosat.org. 

## Contents

* `run_20CRv3_t2m_one_year.sh` - SLURM submissoin shell script that performs a single run on one year (hardcoded from inspection of the public URL link location filename)
* `run_20CRv3_t2m_all_years.sh` - SLURM submission shell script that calls the python script generator code: run_20CRv3_t2m_all_years.py 
* `run_20CRv3_t2m_all_years.py` - python script generator code that loops over all years (hardcoded from inspection of the public URL link location filename) and spawns SLURM submission shell scripts for each year
* `glosat-directory-structure.txt` - Notes from the Confluence WIKI by Jonathan Wynn providing guidance on file management in the GloSAT GWS
* `20CRv3-jasmin.txt` - List of already available 20CRv3 in the GWS
* `20CRv3-t2m-filelist.txt` - List of the 20CRv3 2m Temperature yearly files at source

The first step is to clone the latest jasmin-20CRv3 code into your /home/users/ folder on JASMIN and step into the installed Github directory: 

    $ git clone https://github.com/patternizer/jasmin-20CRv3.git
    $ cd jasmin-20CRv3

### Using Standard Python

The code should run on CentOS-7 systems with the SLURM batch process scheduler and the JasPy module.

Run with:

    $ sbatch run_20CRv3_t2m_one_year.sh (OR)
    $ sbatch run_20CRv3_t2m_all_years.sh 
    
## License

The code is distributed under terms and conditions of the [Open Government License](http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/).

## Contact information

* [Michael Taylor](michael.a.taylor@uea.ac.uk)

