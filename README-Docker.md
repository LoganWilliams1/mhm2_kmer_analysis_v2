Docker support is limited to a single machine and CPUs only (no GPUs).

A MetaHipMer image can be built directly from the bitbucket repository without needing to clone it as follows:

`docker build -t local/mhm2:latest https://bitbucket.org/berkeleylab/mhm2.git#master`

If doing a rebuild, to enure that the latest version is built, use the `--no-cache` option to docker.

To run the docker image, you must set the `--shm-size` parameter; usually 1g to 2g per process will suffice. In addition, you will need to set the mount volume for your inputs and outputs. For example, for a machine with 32 processors:

`docker run --shm-size=32g --volume $datadir:/scratch local/mhm2:latest mhm2.py -r /scratch/arctic_sample_0.fq`

Where `$datadir` is the directory on the host containing the input files and where the output will be written to.
