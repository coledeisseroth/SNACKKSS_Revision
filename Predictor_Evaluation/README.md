#Go into the "import" directory, and gather the needed files specified in BUILD.sh

#Load the Docker image for SNACKKSS_Eval, then run:
docker run -v $(pwd)/:/app/ snackkss-eval-pipeline

