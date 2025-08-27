IMAGE=$1
docker container run -it -p 3838:3838 -v $(pwd)/data:/data ${IMAGE} /bin/bash
