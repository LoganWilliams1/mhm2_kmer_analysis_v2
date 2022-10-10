Docker support is (very) experimental and will ONLY work on a single machine

To build docker, from git:
git submodule init
git submodule update
docker build -f Docker/Dockerfile-build .

To run in docker, you must set the --memory and --shm-size parameters to your specific machine and mount volumes for your inputs and outputs

Examples:

docker run --memory 50g --shm-size=32g --volume /scratch:/scratch robegan21/mhm2:latest upcxx-run -n 24 -shared-heap=1024m /usr/local/bin/mhm2 -r /scratch/arctic_sample_0.fq -o /scratch/test-docker

docker run --memory 50g --shm-size=32g --volume /scratch:/scratch robegan21/mhm2:latest mhm2.py -r /scratch/arctic_sample_0.fq -o /scratch/test-docker
