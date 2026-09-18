# Creating a docker container to run PSOPT

1. Download [Dockerfile](https://github.com/PSOPT/psopt/blob/master/Dockerfile) from the PSOPT distribution, and place it in a folder. This Dockerfile uses [archlinux](https://hub.docker.com/_/archlinux/) as the base. This file clones the latest source code for PSOPT available from GitHub. If you have created your own version (for instance, to include your own examples or cases), you can modify the Dockerfile to copy your own source tree.

2. In your terminal, cd to the same folder where the Dockerfile is. The command to build the docker container (including PSOPT) is as follows: 
```
docker build -t psopt-archlinux:latest .
```
The above command reuses a previous container with the same name, if it exists. If you want to rebuild the whole container use the following command:
```
docker build --no-cache -t psopt-archlinux:latest .
```
3. Issue the following command to run the docker container interactively:
```
docker run -it psopt-archlinux:latest 
```
This will land you in the main 'psopt' folder. From there cd to 'build/examples' to run particular examples, etc.

4. Alternatively, you can use the following command to run the docker container interactively with a data connection to the host 

```
docker run -it --rm -v "$HOME/data:/data" psopt-archlinux:latest 
```
Here, the shared folder is "$HOME/data" as seen from the host, and "/data" as seen from the container.

From within the container, cd to 'build/examples' to run particular examples, etc.
Any output files must be manually copied to the folder /data from within the container. The copied files (e.g. PDFs or .txt files) appear within the corresponding directory of the host ($HOME/data). The host can send files to the container via the same folder.