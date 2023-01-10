# HIV-vaccine-candidate-efficacy
This code accompanies the paper entitled: HIV vaccine candidate efficacy mediated by cyclic AMP-dependent efferocytosis and V2-specific ADCC

To reproduce these results, follow these steps:

1.  Clone this GitHub repo (i.e. the page you are on):
    * ```git clone https://github.com/NIDAP-Community/HIV-vaccine-candidate-efficacy```

2.  The 4 input files should be in the rds_output folder

3.  Install docker and build the docker container:
    * Navigate to the cloned repository. 
    * Move to the ./Docker_file/ directory of this repo

4.  Run the container
    * ```docker build --tag hiv-vaccine-candidate-efficacy```

5. Run the conainer by mounting the ./src/ directory of the repo to /tmp/ in the container and run the following:
    * ```cd /tmp```
    * ```bash run_pipeline.sh```
