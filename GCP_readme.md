## Provisioning a Virtual Machine in Google Cloud
You can create a virtual machine from the existing template: [nemascan-test-vm](https://console.cloud.google.com/compute/instanceTemplates/list?project=andersen-lab). 

Click on **Actions** -> **Create VM**

Scroll to the bottom of the page and click **'Create'**

Once the VM has been created, it should appear in the [list](https://console.cloud.google.com/compute/instances?project=andersen-lab). 

**Remember to Stop the VM and Delete it when you are done to avoid excess charges!!!**

Click on 'SSH' to connect.

Switch to the root user and install some required packages:
```
sudo su
apt-get install docker.io git nano
```


## Using NemaScan Nextflow container
Provision a virtual machine from the existing template and connect with SSH, then switch to the root user and install the required packages:
```
sudo su
apt-get install docker.io git nano
```

To start a terminal session inside of the container, use:
```
docker run -i -t andersenlab/nemascan-nxf /bin/bash
```

Once you are done, exit the container with Ctrl-C or use:
```
exit
```

**Remember to Stop the VM and Delete it when you are done to avoid excess charges!!!**


## Running NemaScan Nextflow container (this is what CeNDR does)
To reproduce the details of a pipeline exactly, you can use the containerized version of the tool set to a specific release version.
Provision a virtual machine from the existing template and connect with SSH, then switch to the root user and install the required packages:
```
sudo su
apt-get install docker.io git nano
```

Use 'nano' to create a *.env file (see: example.env) with variables pointing to your Google Storage locations:
```
nano test.env
```

test.env:
```
TRAIT_FILE="gs://elegansvariation.org/reports/nemascan/abcd123/data.tsv"
OUTPUT_DIR="gs://elegansvariation.org/reports/nemascan/abcd123/results"
WORK_DIR="gs://nf-pipelines/workdir/abcd123"
VCF_VERSION="20210121"
```

```
docker run -i -t \
  --env-file test.env \
  andersenlab/nemascan-nxf:v0.01 \
  nemascan-nxf.sh
```

You can also pass them in as part of the command:
```
docker run -i -t \
  -e TRAIT_FILE="gs://elegansvariation.org/reports/nemascan/abcd123/data.tsv" \
  -e OUTPUT_DIR="gs://elegansvariation.org/reports/nemascan/abcd123/results" \
  -e WORK_DIR="gs://nf-pipelines/workdir/abcd123" \
  -e VCF_VERSION="20210121" \
  andersenlab/nemascan-nxf \
  nemascan-nxf.sh
```

**Remember to Stop the VM and Delete it when you are done to avoid excess charges!!!**



## Testing a new version of the container
First you will have to build and test the container in Google Cloud. Provision a virtual machine from the existing template and connect with SSH, then switch to the root user and install the required packages:
```
sudo su
apt-get install docker.io git nano
```

Clone the repository and check out the branch that contains the version you want to test (in this example, the branch is named 'gcp-nextflow-container'):
```
git clone https://github.com/AndersenLab/NemaScan.git
cd NemaScan
git checkout remotes/origin/gcp-nextflow-container
git pull
```

If you aren't sure of the branch, you can list all available branches with:
```
git branch -a
```

Now you are ready to build the container for testing:
```
docker build -t "andersenlab/nemascan-nxf" . 
```

If the container is built successfully you should be able to see the details with:
```
docker image list
```
example:
```
REPOSITORY                     TAG       IMAGE ID       CREATED          SIZE
andersenlab/nemascan-nxf       latest    bb4f296feec8   26 seconds ago   1.88GB
```


Now you can begin testing. To start the container and open a terminal prompt, substitute your container version's IMAGE ID value in the command below:
```
docker run -i -t bb4f296feec8 /bin/bash
```

Configure the pipeline options by setting the environment variables described in example.env. Substitute the path for your own test data:
```
export TRAIT_FILE="gs://elegansvariation.org/reports/nemascan/abcd123/data.tsv"
export OUTPUT_DIR="gs://elegansvariation.org/reports/nemascan/abcd123/results"
export WORK_DIR="gs://nf-pipelines/workdir/abcd123"
export VCF_VERSION="20210121"
```

Run the pipeline:
```
./nemascan-nxf.sh
```

Because the pipeline takes so long to run, it is possible that your SSH session may time out and disconnect during the test.
Reconnect to the VM with SSH and then list the running containers with:
```
docker container list
```
You should see an output similar to this:
```
CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS              PORTS               NAMES
d2f3dbb5e136        bb4f296feec8        "/bin/bash"         6 minutes ago       Up 6 minutes                            jolly_cartwright
```

You can re-attach to the container by substituting your own CONTAINER ID in the command below:
```
docker attach d2f3dbb5e136
```
Warning: Nextflow only prints status updates sporadically, so you may not see any output for some time after attaching to the container.

After the pipeline has completed, use Ctrl-C to exit the container, or type:
```
exit
```


## Publishing a new version of the container
Once you have validated the pipeline against test data, you can publish a release to docker hub.
Switch to the root user and log in with your credentials to docker hub:

```
sudo su
docker login
```

List the containers, and select the IMAGE ID of the container you just tested:
```
docker image list
```

example:
```
REPOSITORY                 TAG                 IMAGE ID            CREATED             SIZE
andersenlab/nemascan-nxf   latest              bb4f296feec8        4 hours ago         1.88GB
google/cloud-sdk           slim                d6d0a7854ac3        28 hours ago        1.16GB
```

Substitute your IMAGE ID, then tag the container with the docker hub repository, container name, and version number. You can also omit the version number to make the selected version the default (or 'latest')
```
docker tag bb4f296feec8 andersenlab/nemascan-nxf
```
```
docker tag bb4f296feec8 andersenlab/nemascan-nxf:v0.01
```

Publish the container to docker hub with a version tag:
```
docker push andersenlab/nemascan-nxf:v0.01
```
You can also publish the container to docker hub without a version number if you want it to be the default ('latest') version:
```
docker push andersenlab/nemascan-nxf
```

**Remember to Stop the VM and Delete it when you are done to avoid excess charges!!!**
