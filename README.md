# HIV-vaccine-candidate-efficacy
This code accompanies the paper: Bissa M, Kim S, Galli V, Fourati S, Sarkis S, Arakelyan A, de Castro IS, Rahman MA, Fujiwara S, Vaccari M, Tomalka JA, Stamos JD, Schifanella L, Gorini G, Moles R, Gutowska A, Ferrari G, Lobanov A, Montefiori DC, Nelson GW, Cam MC, Chakhtoura M, Haddad EK, Doster MN, McKinnon K, Brown S, Venzon DJ, Choo-Wosoba H, Breed MW, Killoran KE, Kramer J, Margolis L, Sekaly RP, Hager GL, Franchini G. HIV vaccine candidate efficacy in female macaques mediated by cAMP-dependent efferocytosis and V2-specific ADCC. Nat Commun. 2023 Feb 2;14(1):575. doi: 10.1038/s41467-023-36109-8. [Pubmed Link](https://pubmed.ncbi.nlm.nih.gov/36732510/)

To reproduce these results, follow these steps:

1.  Clone this GitHub repo (i.e. the page you are on):
    * ```git clone https://github.com/NIDAP-Community/HIV-vaccine-candidate-efficacy.git```

2.  The 4 input files should be in the ./src/rds_output folder

3.  Install docker and build the docker container:
    * Navigate to the cloned repository directory. 
    * Move to the ./Docker_file/ directory of this repo

4.  Build the container:
    * ```docker build --tag hiv-vaccine-candidate-efficacy .```

5.  Navidate to the cloned repository directory, Run the conainer by mounting the ./src/ directory of the repo to /tmp/ in the container:
    * ```docker run -ti -v $(pwd)/src:/tmp hiv-vaccine-candidate-efficacy```

6.  Activate the conda environment in the docker container:
    * ```conda activate legacy-default```

5.  Run the following code
    * ```cd /tmp```
    * ```bash run_pipeline.sh```
