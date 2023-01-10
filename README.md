# HIV-vaccine-candidate-efficacy
This code accompanies the paper entitled: HIV vaccine candidate efficacy mediated by cyclic AMP-dependent efferocytosis and V2-specific ADCC

To reproduce these results, follow these steps:

Clone this GitHub repo (i.e. the page you are on):
git clone https://github.com/NIDAP-Community/HIV-vaccine-candidate-efficacy
The 4 input files should be in the rds_output folder

Install docker and build the docker container:
Navigate to the cloned repository. Move to the ./Docker_file/ directory of this repo
Run: docker build --tag hiv-vaccine-candidate-efficacy
Run the conainer by mounting the ./src/ directory of the repo to /tmp/ in the container and run the following:
docker run -ti -v $(pwd)/src:/tmp minnar_el_al_2022
cd /tmp
bash run_pipeline.sh
