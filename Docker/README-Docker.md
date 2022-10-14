Docker support is (very) experimental and will ONLY work on a single machine.

To build docker, from git:
git submodule init
git submodule update
docker build -t local/mhm2:latest -f Docker/Dockerfile-build .

To run in docker, you must set the --memory and --shm-size parameters to your specific machine and mount volumes for your inputs and outputs

Example:

docker run --memory 50g --shm-size=32g --volume $datadir:/scratch local/mhm2:latest mhm2.py -r /scratch/arctic_sample_0.fq -o /scratch/test-docker

Where $datadir is the directory on the host containing the input files and where the output will be written to.
